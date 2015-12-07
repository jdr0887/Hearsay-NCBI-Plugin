package org.renci.hearsay.commands.ncbi;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.Scanner;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import org.apache.commons.collections4.CollectionUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.math.NumberUtils;
import org.renci.gbff.GBFFFilter;
import org.renci.gbff.GBFFManager;
import org.renci.gbff.filter.GBFFAndFilter;
import org.renci.gbff.filter.GBFFFeatureSourceOrganismNameFilter;
import org.renci.gbff.filter.GBFFFeatureTypeNameFilter;
import org.renci.gbff.filter.GBFFSequenceAccessionPrefixFilter;
import org.renci.gbff.filter.GBFFSourceOrganismNameFilter;
import org.renci.gbff.model.Feature;
import org.renci.gbff.model.Sequence;
import org.renci.hearsay.commands.ncbi.util.FTPUtil;
import org.renci.hearsay.dao.HearsayDAOBeanService;
import org.renci.hearsay.dao.model.Alignment;
import org.renci.hearsay.dao.model.Chromosome;
import org.renci.hearsay.dao.model.Gene;
import org.renci.hearsay.dao.model.Identifier;
import org.renci.hearsay.dao.model.Location;
import org.renci.hearsay.dao.model.ReferenceSequence;
import org.renci.hearsay.dao.model.Region;
import org.renci.hearsay.dao.model.RegionType;
import org.renci.hearsay.dao.model.StrandType;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class PullAlignmentsRunnable implements Runnable {

    private static final Logger logger = LoggerFactory.getLogger(PullAlignmentsRunnable.class);

    private static final Pattern featureLocationPattern = Pattern.compile("^(join|order)\\((.+)\\)$");

    private static final List<String> inclusionPatterns = Arrays.asList(new String[] { "misc_feature", "polyA_signal", "polyA_site",
            "transit_peptide", "mat_peptide", "sig_peptide", "unsure", "stem_loop", "protein_bind", "repeat_region", "prim_transcript",
            "proprotein", "LTR", "TATA_signal", "primer_bind", "terminator", "misc_difference", "misc_binding", "RBS", "misc_signal",
            "J_segment", "C_region", "conflict", "promoter", "ncRNA", "modified_base" });

    private HearsayDAOBeanService hearsayDAOBeanService;

    public PullAlignmentsRunnable(HearsayDAOBeanService hearsayDAOBeanService) {
        super();
        this.hearsayDAOBeanService = hearsayDAOBeanService;
    }

    @Override
    public void run() {
        logger.debug("ENTERING run()");

        try {

            LinkedList<String> alignmentsLines = new LinkedList<String>();

            File alignmentsFile = FTPUtil.ncbiDownload("/refseq/H_sapiens/alignments", "GCF_000001405.28_knownrefseq_alignments.gff3");
            try (FileInputStream fis = new FileInputStream(alignmentsFile);
                    InputStreamReader isr = new InputStreamReader(fis);
                    BufferedReader br = new BufferedReader(isr)) {
                String line;
                while ((line = br.readLine()) != null) {
                    alignmentsLines.add(line);
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
            alignmentsFile.delete();

            // this will take a while
            GBFFManager gbffMgr = GBFFManager.getInstance(1, true);

            List<GBFFFilter> filters = Arrays
                    .asList(new GBFFFilter[] { new GBFFSequenceAccessionPrefixFilter(Arrays.asList(new String[] { "NM_", "NR_" })),
                            new GBFFSourceOrganismNameFilter("Homo sapiens"), new GBFFFeatureSourceOrganismNameFilter("Homo sapiens"),
                            new GBFFFeatureTypeNameFilter("CDS"), new GBFFFeatureTypeNameFilter("source") });

            GBFFAndFilter gbffFilter = new GBFFAndFilter(filters);

            List<File> fileList = FTPUtil.ncbiDownloadBySuffix("/refseq/H_sapiens/mRNA_Prot", "rna.gbff.gz");

            ExecutorService es = Executors.newFixedThreadPool(8);
            for (File f : fileList) {
                List<Sequence> sequenceList = gbffMgr.deserialize(gbffFilter, f);
                logger.debug("sequenceList.size(): {}", sequenceList.size());
                if (CollectionUtils.isNotEmpty(sequenceList)) {

                    for (Sequence sequence : sequenceList) {

                        String refSeqVersionedAccession = sequence.getVersion().trim().contains(" ")
                                ? sequence.getVersion().substring(0, sequence.getVersion().indexOf(" ")) : sequence.getVersion();

                        List<Identifier> rnaNucleotideAccessionIdentifierList = hearsayDAOBeanService.getIdentifierDAO()
                                .findByExample(new Identifier("www.ncbi.nlm.nih.gov/nuccore", refSeqVersionedAccession));

                        String proteinAccession = null;
                        Feature firstCDSFeature = null;
                        for (Feature feature : sequence.getFeatures()) {
                            if (!"CDS".equals(feature.getType())) {
                                continue;
                            }
                            firstCDSFeature = feature;
                            break;
                        }
                        proteinAccession = firstCDSFeature.getQualifiers().getProperty("protein_id").replace("\"", "");

                        List<Identifier> proteinAccessionIdentifierList = hearsayDAOBeanService.getIdentifierDAO()
                                .findByExample(new Identifier("www.ncbi.nlm.nih.gov/protein", proteinAccession));

                        List<Long> identifierIdList = new ArrayList<Long>();
                        for (Identifier identifier : rnaNucleotideAccessionIdentifierList) {
                            identifierIdList.add(identifier.getId());
                        }
                        for (Identifier identifier : proteinAccessionIdentifierList) {
                            identifierIdList.add(identifier.getId());
                        }

                        List<ReferenceSequence> potentialRefSeqs = hearsayDAOBeanService.getReferenceSequenceDAO()
                                .findByIdentifiers(identifierIdList);

                        if (CollectionUtils.isEmpty(potentialRefSeqs)) {
                            logger.warn("Could not find ReferenceSequence: refSeqVersionedAccession = {}, proteinAccession = {}",
                                    refSeqVersionedAccession, proteinAccession);
                            continue;
                        }

                        logger.info("Using ReferenceSequence: refSeqVersionedAccession = {}, proteinAccession = {}",
                                refSeqVersionedAccession, proteinAccession);

                        // persistAlignmentsExecutorService.submit(new PersistAlignmentsFromUCSCRunnable(hearsayDAOBean,
                        // sequence));
                        es.submit(new PersistAlignmentsFromNCBIRunnable(sequence, potentialRefSeqs, firstCDSFeature, alignmentsLines));
                        es.submit(new PersistFeaturesRunnable(sequence, potentialRefSeqs));
                    }

                }
                // f.delete();
            }
            es.shutdown();
            es.awaitTermination(4L, TimeUnit.HOURS);

        } catch (Exception e) {
            logger.error(e.getMessage(), e);
        }

    }

    class PersistAlignmentsFromNCBIRunnable implements Runnable {

        private final Sequence sequence;

        private final List<ReferenceSequence> referenceSequences;

        private final Feature firstCDSFeature;

        private final LinkedList<String> alignmentsLines;

        public PersistAlignmentsFromNCBIRunnable(final Sequence sequence, final List<ReferenceSequence> referenceSequences,
                final Feature firstCDSFeature, final LinkedList<String> alignmentsLines) {
            super();
            this.sequence = sequence;
            this.referenceSequences = referenceSequences;
            this.firstCDSFeature = firstCDSFeature;
            this.alignmentsLines = alignmentsLines;
        }

        @Override
        public void run() {
            logger.debug("ENTERING run()");

            try {
                String refSeqVersionedAccession = sequence.getVersion().trim().contains(" ")
                        ? sequence.getVersion().substring(0, sequence.getVersion().indexOf(" ")) : sequence.getVersion();

                final StrandType strandType = referenceSequences.get(0).getStrandType();

                // add protein info to alignment
                String firstCDSFeatureLocation = firstCDSFeature.getLocation();
                Integer firstCDSFeatureLocationStart = null;
                Integer firstCDSFeatureLocationStop = null;
                if (firstCDSFeatureLocation.contains("join")) {
                    Matcher m = featureLocationPattern.matcher(firstCDSFeatureLocation);
                    m.find();
                    String joinContent = m.group(2);
                    List<Integer> positions = new ArrayList<Integer>();
                    String[] ranges = joinContent.split(",");
                    for (String r : ranges) {
                        String[] split = r.split("\\.\\.");
                        positions.add(Integer.valueOf(split[0]));
                        positions.add(Integer.valueOf(split[1]));
                    }

                    Collections.sort(positions);
                    firstCDSFeatureLocationStart = positions.get(0);
                    firstCDSFeatureLocationStop = positions.get(positions.size() - 1);

                } else {
                    String[] split = firstCDSFeatureLocation.split("\\.\\.");
                    if (NumberUtils.isNumber(split[0]) && NumberUtils.isNumber(split[1])) {
                        firstCDSFeatureLocationStart = Integer.valueOf(split[0]);
                        firstCDSFeatureLocationStop = Integer.valueOf(split[1]);
                    }
                }

                int exonCount = 0;
                if (CollectionUtils.isNotEmpty(sequence.getFeatures())) {
                    for (Feature feature : sequence.getFeatures()) {
                        if (!"exon".equals(feature.getType())) {
                            continue;
                        }
                        exonCount++;
                    }
                }
                if (exonCount == 0) {
                    logger.warn("empty features");
                    return;
                }

                // add alignments
                Alignment alignment = new Alignment();
                alignment.setId(hearsayDAOBeanService.getAlignmentDAO().save(alignment));
                alignment.getReferenceSequences().addAll(referenceSequences);

                Location proteinLocation = new Location(firstCDSFeatureLocationStart, firstCDSFeatureLocationStop);
                proteinLocation.setId(hearsayDAOBeanService.getLocationDAO().save(proteinLocation));
                alignment.setProteinLocation(proteinLocation);

                // add exons to alignment
                for (Feature feature : sequence.getFeatures()) {

                    if (!"exon".equals(feature.getType())) {
                        continue;
                    }

                    Region region = new Region(RegionType.EXON);
                    region.setAlignment(alignment);
                    String range = feature.getLocation();
                    String start = null;
                    String stop = null;
                    String[] split = range.split("\\.\\.");

                    if (strandType.equals(StrandType.MINUS)) {
                        stop = split[0];
                        start = split[1];
                    } else {
                        stop = split[1];
                        start = split[0];
                    }

                    if (NumberUtils.isNumber(start) && NumberUtils.isNumber(stop)) {
                        Location transcriptLocation = new Location(Integer.valueOf(start), Integer.valueOf(stop));
                        transcriptLocation.setId(hearsayDAOBeanService.getLocationDAO().save(transcriptLocation));
                        region.setTranscriptLocation(transcriptLocation);
                    }

                    region.setId(hearsayDAOBeanService.getRegionDAO().save(region));
                    alignment.getRegions().add(region);
                }

                hearsayDAOBeanService.getAlignmentDAO().save(alignment);

                List<Region> regions = alignment.getRegions();
                if (CollectionUtils.isNotEmpty(regions)) {
                    Collections.sort(regions, new Comparator<Region>() {
                        @Override
                        public int compare(Region r1, Region r2) {
                            if (strandType.equals(StrandType.MINUS)) {
                                return r2.getTranscriptLocation().getStart().compareTo(r1.getTranscriptLocation().getStart());
                            } else {
                                return r1.getTranscriptLocation().getStart().compareTo(r2.getTranscriptLocation().getStart());
                            }
                        }
                    });

                    String refSeqGenomicAccession = null;
                    for (Identifier identifier : referenceSequences.get(0).getIdentifiers()) {
                        if (identifier.getSystem().equals("www.ncbi.nlm.nih.gov/genome")) {
                            refSeqGenomicAccession = identifier.getValue();
                            break;
                        }
                    }

                    logger.debug("refSeqGenomicAccession = {}", refSeqGenomicAccession);
                    for (String line : alignmentsLines) {

                        if (line.startsWith("#")) {
                            continue;
                        }

                        // NT_187633.1 RefSeq cDNA_match 278295 278486 192 - . ID=aln45579;Target=NM_000853.3 1 192
                        // +;gap_count=0;identity=1;idty=1;num_ident=1093;num_mismatch=0;pct_coverage=100;pct_identity_gap=100;pct_identity_ungap=100;score=192
                        String[] columns = line.split("\t");
                        String genomicAccession = columns[0];
                        String genomicStart = columns[3];
                        String genomicStop = columns[4];
                        String strand = columns[6];
                        String attributes = columns[8];

                        if (!genomicAccession.equals(refSeqGenomicAccession)) {
                            continue;
                        }

                        if (attributes.contains(refSeqVersionedAccession)) {

                            String[] attributeSplit = attributes.split(";");
                            String[] targetSplit = attributeSplit[1].split(" ");
                            Integer start = Integer.valueOf(targetSplit[1]);
                            Integer stop = Integer.valueOf(targetSplit[2]);

                            for (Region region : regions) {
                                if (strandType.equals(StrandType.MINUS)) {
                                    if (region.getTranscriptLocation().getStart().equals(stop)
                                            && region.getTranscriptLocation().getStop().equals(start)) {
                                        Location genomicLocation = new Location(Integer.valueOf(genomicStart),
                                                Integer.valueOf(genomicStop));
                                        genomicLocation.setId(hearsayDAOBeanService.getLocationDAO().save(genomicLocation));
                                        region.setRegionLocation(genomicLocation);
                                        hearsayDAOBeanService.getRegionDAO().save(region);
                                    }
                                } else {
                                    if (region.getTranscriptLocation().getStart().equals(start)
                                            && region.getTranscriptLocation().getStop().equals(stop)) {
                                        Location genomicLocation = new Location(Integer.valueOf(genomicStart),
                                                Integer.valueOf(genomicStop));
                                        genomicLocation.setId(hearsayDAOBeanService.getLocationDAO().save(genomicLocation));
                                        region.setRegionLocation(genomicLocation);
                                        hearsayDAOBeanService.getRegionDAO().save(region);
                                    }
                                }
                            }

                        }

                    }
                }

            } catch (Exception e) {
                logger.error(e.getMessage(), e);
                e.printStackTrace();
            }

        }

    }

    class PersistAlignmentsFromUCSCRunnable implements Runnable {

        private final Sequence sequence;

        public PersistAlignmentsFromUCSCRunnable(final Sequence sequence) {
            super();
            this.sequence = sequence;
        }

        @Override
        public void run() {
            logger.debug("ENTERING run()");

            LinkedList<String> refGeneLines = new LinkedList<String>();

            // 5MB gzipped
            File refGeneFile = FTPUtil.ucscDownload("/goldenPath/hg38/database", "refGene.txt.gz");
            try (FileInputStream fis = new FileInputStream(refGeneFile);
                    GZIPInputStream gis = new GZIPInputStream(fis);
                    InputStreamReader isr = new InputStreamReader(gis);
                    BufferedReader br = new BufferedReader(isr)) {
                String line;
                while ((line = br.readLine()) != null) {
                    refGeneLines.add(line);
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
            // refGeneFile.delete();

            try {
                String refseqAccession = sequence.getAccession().contains(" ")
                        ? sequence.getAccession().substring(0, sequence.getAccession().indexOf(" ")) : sequence.getAccession();

                String refSeqVersionedAccession = sequence.getVersion().trim().contains(" ")
                        ? sequence.getVersion().substring(0, sequence.getVersion().indexOf(" ")) : sequence.getVersion();

                List<Identifier> rnaNucleotideAccessionIdentifierList = hearsayDAOBeanService.getIdentifierDAO()
                        .findByExample(new Identifier("www.ncbi.nlm.nih.gov/nuccore", refSeqVersionedAccession));

                String proteinAccession = null;
                Feature firstCDSFeature = null;
                for (Feature feature : sequence.getFeatures()) {
                    if (!"CDS".equals(feature.getType())) {
                        continue;
                    }
                    firstCDSFeature = feature;
                    break;
                }

                proteinAccession = firstCDSFeature.getQualifiers().getProperty("protein_id").replace("\"", "");

                List<Identifier> proteinAccessionIdentifierList = hearsayDAOBeanService.getIdentifierDAO()
                        .findByExample(new Identifier("www.ncbi.nlm.nih.gov/protein", proteinAccession));

                List<Long> identifierIdList = new ArrayList<Long>();
                for (Identifier identifier : rnaNucleotideAccessionIdentifierList) {
                    identifierIdList.add(identifier.getId());
                }
                for (Identifier identifier : proteinAccessionIdentifierList) {
                    identifierIdList.add(identifier.getId());
                }

                List<ReferenceSequence> potentialRefSeqs = hearsayDAOBeanService.getReferenceSequenceDAO()
                        .findByIdentifiers(identifierIdList);

                if (CollectionUtils.isEmpty(potentialRefSeqs)) {
                    logger.warn("Could not find ReferenceSequence: refSeqVersionedAccession = {}, proteinAccession = {}",
                            refSeqVersionedAccession, proteinAccession);
                    return;
                }

                logger.info("Using ReferenceSequence: refSeqVersionedAccession = {}, proteinAccession = {}", refSeqVersionedAccession,
                        proteinAccession);

                final StrandType strandType = potentialRefSeqs.get(0).getStrandType();

                // add alignments
                Alignment alignment = new Alignment();
                alignment.setId(hearsayDAOBeanService.getAlignmentDAO().save(alignment));
                alignment.getReferenceSequences().addAll(potentialRefSeqs);

                // add protein info to alignment
                String firstCDSFeatureLocation = firstCDSFeature.getLocation();
                Integer firstCDSFeatureLocationStart = null;
                Integer firstCDSFeatureLocationStop = null;
                if (firstCDSFeatureLocation.contains("join")) {
                    Matcher m = featureLocationPattern.matcher(firstCDSFeatureLocation);
                    m.find();
                    String joinContent = m.group(2);
                    List<Integer> positions = new ArrayList<Integer>();
                    String[] ranges = joinContent.split(",");
                    for (String r : ranges) {
                        String[] split = r.split("\\.\\.");
                        positions.add(Integer.valueOf(split[0]));
                        positions.add(Integer.valueOf(split[1]));
                    }

                    Collections.sort(positions);
                    firstCDSFeatureLocationStart = positions.get(0);
                    firstCDSFeatureLocationStop = positions.get(positions.size() - 1);

                } else {
                    String[] split = firstCDSFeatureLocation.split("\\.\\.");
                    if (NumberUtils.isNumber(split[0]) && NumberUtils.isNumber(split[1])) {
                        firstCDSFeatureLocationStart = Integer.valueOf(split[0]);
                        firstCDSFeatureLocationStop = Integer.valueOf(split[1]);
                    }
                }

                Location proteinLocation = new Location(firstCDSFeatureLocationStart, firstCDSFeatureLocationStop);
                proteinLocation.setId(hearsayDAOBeanService.getLocationDAO().save(proteinLocation));
                alignment.setProteinLocation(proteinLocation);

                hearsayDAOBeanService.getAlignmentDAO().save(alignment);

                // add exons to alignment
                List<Region> regions = new ArrayList<Region>();
                for (Feature feature : sequence.getFeatures()) {
                    if (!"exon".equals(feature.getType())) {
                        continue;
                    }

                    Region region = new Region(RegionType.EXON);
                    String range = feature.getLocation();
                    String start = null;
                    String stop = null;
                    String[] split = range.split("\\.\\.");

                    if (strandType.equals(StrandType.MINUS)) {
                        stop = split[0];
                        start = split[1];
                    } else {
                        stop = split[1];
                        start = split[0];
                    }

                    if (NumberUtils.isNumber(start) && NumberUtils.isNumber(stop)) {
                        Location transcriptLocation = new Location(Integer.valueOf(start), Integer.valueOf(stop));
                        transcriptLocation.setId(hearsayDAOBeanService.getLocationDAO().save(transcriptLocation));
                        logger.debug(transcriptLocation.toString());
                        region.setTranscriptLocation(transcriptLocation);
                    }
                    regions.add(region);
                }

                regions.sort((a, b) -> {
                    if (strandType.equals(StrandType.MINUS)) {
                        return b.getTranscriptLocation().getStart().compareTo(a.getTranscriptLocation().getStart());
                    } else {
                        return a.getTranscriptLocation().getStart().compareTo(b.getTranscriptLocation().getStart());
                    }
                });

                Gene refseqGene = potentialRefSeqs.get(0).getGene();
                List<String> chromosomeList = new ArrayList<String>();
                for (Chromosome chr : refseqGene.getChromosomes()) {
                    String chromosomeName = chr.getName();
                    chromosomeList.add(chromosomeName);
                }

                List<Location> genomicLocationList = new ArrayList<Location>();
                for (String line : refGeneLines) {
                    try (Scanner scanner = new Scanner(line).useDelimiter("\t")) {
                        String bin = scanner.next();
                        String name = scanner.next();
                        if (!name.equals(refseqAccession)) {
                            continue;
                        }
                        String chromosome = scanner.next();

                        if (!chromosomeList.contains(chromosome.replace("chr", ""))) {
                            continue;
                        }

                        String strand = scanner.next();
                        if (!strand.equals(strandType.getValue())) {
                            continue;
                        }
                        String txStart = scanner.next();
                        String txEnd = scanner.next();
                        String cdsStart = scanner.next();
                        String cdsEnd = scanner.next();
                        String exonCount = scanner.next();

                        String exonStarts = scanner.next();
                        String[] exonStartsSplit = exonStarts.split(",");

                        String exonEnds = scanner.next();
                        String[] exonEndsSplit = exonEnds.split(",");

                        for (int i = 0; i < Integer.valueOf(exonCount); i++) {
                            Location genomicLocation = new Location(Integer.valueOf(exonStartsSplit[i]), Integer.valueOf(exonEndsSplit[i]));
                            logger.debug(genomicLocation.toString());
                            genomicLocationList.add(genomicLocation);
                        }
                        scanner.close();
                    }
                }

                if (regions.size() != genomicLocationList.size()) {
                    logger.warn("regions.size() = {}, genomicLocationList.size() = {}", regions.size(), genomicLocationList.size());
                    return;
                }

                genomicLocationList.sort((a, b) -> a.getStart().compareTo(b.getStart()));

                // initial exon regions
                for (int i = 0; i < regions.size(); i++) {
                    Region region = regions.get(i);
                    region.setAlignment(alignment);
                    Location regionLocation = genomicLocationList.get(i);
                    regionLocation.setId(hearsayDAOBeanService.getLocationDAO().save(regionLocation));
                    region.setRegionLocation(regionLocation);
                    region.setId(hearsayDAOBeanService.getRegionDAO().save(region));
                    alignment.getRegions().add(region);
                }

                hearsayDAOBeanService.getAlignmentDAO().save(alignment);

            } catch (Exception e) {
                logger.error(e.getMessage(), e);
                e.printStackTrace();
            }

        }

    }

    class PersistFeaturesRunnable implements Runnable {

        private final Sequence sequence;

        private final List<ReferenceSequence> referenceSequences;

        public PersistFeaturesRunnable(final Sequence sequence, final List<ReferenceSequence> referenceSequences) {
            super();
            this.sequence = sequence;
            this.referenceSequences = referenceSequences;
        }

        @Override
        public void run() {
            logger.debug("ENTERING run()");

            try {

                // add features
                for (Feature feature : sequence.getFeatures()) {
                    if (!inclusionPatterns.contains(feature.getType())) {
                        continue;
                    }
                    logger.debug(feature.toString());
                    String location = feature.getLocation();
                    org.renci.hearsay.dao.model.Feature hearsayFeature = new org.renci.hearsay.dao.model.Feature();
                    hearsayFeature.setType(feature.getType());
                    String note = feature.getQualifiers().getProperty("note");

                    if (StringUtils.isNotEmpty(note)) {
                        hearsayFeature.setNote(note);
                    }

                    hearsayFeature.setId(hearsayDAOBeanService.getFeatureDAO().save(hearsayFeature));
                    hearsayFeature.getReferenceSequences().addAll(referenceSequences);

                    if (NumberUtils.isNumber(location)) {
                        Location l = new Location(Integer.valueOf(location), Integer.valueOf(location));
                        l.setId(hearsayDAOBeanService.getLocationDAO().save(l));
                        hearsayFeature.getLocations().add(l);
                    } else if (location.startsWith("join") || location.startsWith("order")) {

                        Matcher m = featureLocationPattern.matcher(location);
                        m.find();
                        try (Scanner scanner = new Scanner(m.group(2)).useDelimiter(",")) {
                            while (scanner.hasNext()) {
                                String range = scanner.next();
                                String startValue = range.substring(0, range.indexOf(".."));
                                String stopValue = range.substring(range.indexOf("..") + 2, range.length());
                                if (NumberUtils.isNumber(startValue) && NumberUtils.isNumber(stopValue)) {
                                    Location l = new Location(Integer.valueOf(startValue), Integer.valueOf(stopValue));
                                    l.setId(hearsayDAOBeanService.getLocationDAO().save(l));
                                    hearsayFeature.getLocations().add(l);
                                }
                            }
                            scanner.close();
                        }

                    } else if (location.contains("..")) {
                        String startValue = location.substring(0, location.indexOf(".."));
                        String stopValue = location.substring(location.indexOf("..") + 2, location.length());
                        if (NumberUtils.isNumber(startValue) && NumberUtils.isNumber(stopValue)) {
                            Location l = new Location(Integer.valueOf(startValue), Integer.valueOf(stopValue));
                            l.setId(hearsayDAOBeanService.getLocationDAO().save(l));
                            hearsayFeature.getLocations().add(l);
                        }
                    }

                    hearsayDAOBeanService.getFeatureDAO().save(hearsayFeature);
                    logger.debug(hearsayFeature.toString());

                }
            } catch (Exception e) {
                logger.error(e.getMessage(), e);
                e.printStackTrace();
            }

        }

    }

}
