package org.renci.hearsay.commands.ncbi;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Properties;
import java.util.Scanner;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

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
import org.renci.gff3.GFF3Manager;
import org.renci.gff3.filters.AttributeValueFilter;
import org.renci.gff3.model.GFF3Record;
import org.renci.hearsay.commands.ncbi.util.FTPUtil;
import org.renci.hearsay.dao.HearsayDAOBeanService;
import org.renci.hearsay.dao.HearsayDAOException;
import org.renci.hearsay.dao.model.Alignment;
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

            File alignmentsFile = FTPUtil.ncbiDownload("/refseq/H_sapiens/alignments", "GCF_000001405.28_knownrefseq_alignments.gff3");
            GFF3Manager gff3Mgr = GFF3Manager.getInstance(alignmentsFile);

            // this will take a while
            GBFFManager gbffMgr = GBFFManager.getInstance(1, true);

            List<GBFFFilter> filters = Arrays
                    .asList(new GBFFFilter[] { new GBFFSequenceAccessionPrefixFilter(Arrays.asList(new String[] { "NM_", "NR_" })),
                            new GBFFSourceOrganismNameFilter("Homo sapiens"), new GBFFFeatureSourceOrganismNameFilter("Homo sapiens"),
                            new GBFFFeatureTypeNameFilter("CDS"), new GBFFFeatureTypeNameFilter("source") });

            GBFFAndFilter gbffFilter = new GBFFAndFilter(filters);

            List<File> fileList = FTPUtil.ncbiDownloadBySuffix("/refseq/H_sapiens/mRNA_Prot", "rna.gbff.gz");

            List<Sequence> allSequences = new ArrayList<Sequence>();

            for (File f : fileList) {
                List<Sequence> sequenceList = gbffMgr.deserialize(gbffFilter, f);
                if (CollectionUtils.isNotEmpty(sequenceList)) {
                    logger.debug("sequenceList.size(): {}", sequenceList.size());
                    allSequences.addAll(sequenceList);
                }
            }

            ExecutorService es = Executors.newFixedThreadPool(4);

            if (CollectionUtils.isNotEmpty(allSequences)) {

                for (Sequence sequence : allSequences) {

                    List<Identifier> identifierList = new ArrayList<Identifier>();

                    // rna nucleotide accession
                    String refSeqVersionedAccession = sequence.getVersion().trim().contains(" ")
                            ? sequence.getVersion().substring(0, sequence.getVersion().indexOf(" ")) : sequence.getVersion();
                    List<Identifier> rnaNucleotideAccessionIdentifierList = hearsayDAOBeanService.getIdentifierDAO()
                            .findByExample(new Identifier("www.ncbi.nlm.nih.gov/nuccore", refSeqVersionedAccession));
                    if (CollectionUtils.isNotEmpty(rnaNucleotideAccessionIdentifierList)) {
                        identifierList.add(rnaNucleotideAccessionIdentifierList.get(0));
                    }

                    // protein accession
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
                    if (CollectionUtils.isNotEmpty(proteinAccessionIdentifierList)) {
                        identifierList.add(rnaNucleotideAccessionIdentifierList.get(0));
                    }

                    identifierList.forEach(a -> logger.info(a.toString()));

                    List<ReferenceSequence> potentialRefSeqs = hearsayDAOBeanService.getReferenceSequenceDAO()
                            .findByIdentifiers(identifierList);

                    if (CollectionUtils.isEmpty(potentialRefSeqs)) {
                        logger.warn("Could not find ReferenceSequence");
                        continue;
                    }

                    ReferenceSequence referenceSequence = potentialRefSeqs.get(0);

                    // submitting
                    es.submit(new Runnable() {

                        @Override
                        public void run() {
                            try {

                                // add features
                                for (Feature feature : sequence.getFeatures()) {
                                    if (!inclusionPatterns.contains(feature.getType())) {
                                        continue;
                                    }
                                    logger.debug(feature.toString());
                                    org.renci.hearsay.dao.model.Feature hearsayFeature = new org.renci.hearsay.dao.model.Feature(
                                            feature.getType());
                                    String note = feature.getQualifiers().getProperty("note");
                                    if (StringUtils.isNotEmpty(note)) {
                                        hearsayFeature.setNote(note);
                                    }
                                    hearsayFeature.setId(hearsayDAOBeanService.getFeatureDAO().save(hearsayFeature));
                                    hearsayFeature.getReferenceSequences().add(referenceSequence);

                                    String location = feature.getLocation();
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
                            }
                        }
                    });

                    String firstCDSFeatureLocation = firstCDSFeature.getLocation();

                    es.submit(new Runnable() {

                        @Override
                        public void run() {

                            try {

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
                                    logger.warn("no exons found");
                                    return;
                                }

                                String refSeqVersionedAccession = sequence.getVersion().trim().contains(" ")
                                        ? sequence.getVersion().substring(0, sequence.getVersion().indexOf(" ")) : sequence.getVersion();

                                StrandType strandType = referenceSequence.getStrandType();

                                // add protein info to alignment
                                Location proteinLocation = null;
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
                                    proteinLocation = new Location(positions.get(0), positions.get(positions.size() - 1));
                                } else {
                                    String[] split = firstCDSFeatureLocation.split("\\.\\.");
                                    if (NumberUtils.isNumber(split[0]) && NumberUtils.isNumber(split[1])) {
                                        proteinLocation = new Location(Integer.valueOf(split[0]), Integer.valueOf(split[1]));
                                    }
                                }

                                // add alignments
                                Alignment alignment = new Alignment();
                                if (proteinLocation != null) {
                                    proteinLocation.setId(hearsayDAOBeanService.getLocationDAO().save(proteinLocation));
                                    alignment.setProteinLocation(proteinLocation);
                                }
                                alignment.getReferenceSequences().add(referenceSequence);

                                // add exons to alignment
                                for (Feature feature : sequence.getFeatures()) {

                                    if (!"exon".equals(feature.getType())) {
                                        continue;
                                    }

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

                                        Region region = new Region(RegionType.EXON);
                                        region.setTranscriptLocation(transcriptLocation);
                                        region.setAlignment(alignment);
                                        region.setId(hearsayDAOBeanService.getRegionDAO().save(region));

                                        alignment.getRegions().add(region);
                                    }
                                }
                                alignment.setId(hearsayDAOBeanService.getAlignmentDAO().save(alignment));

                                List<Region> regions = alignment.getRegions();
                                if (CollectionUtils.isNotEmpty(regions)) {
                                    if (strandType.equals(StrandType.MINUS)) {
                                        regions.sort((a, b) -> b.getTranscriptLocation().getStart()
                                                .compareTo(a.getTranscriptLocation().getStart()));
                                    } else {
                                        regions.sort((a, b) -> a.getTranscriptLocation().getStart()
                                                .compareTo(b.getTranscriptLocation().getStart()));
                                    }

                                    List<GFF3Record> records = gff3Mgr
                                            .deserialize(new AttributeValueFilter("Target", refSeqVersionedAccession));

                                    records.forEach(a -> {

                                        Properties attributes = a.getAttributes();
                                        String targetValue = attributes.get("Target").toString();
                                        String[] targetSplit = targetValue.split(" ");

                                        Integer start = Integer.valueOf(targetSplit[1]);
                                        Integer stop = Integer.valueOf(targetSplit[2]);

                                        try {
                                            for (Region region : regions) {
                                                if (strandType.equals(StrandType.MINUS)) {
                                                    if (region.getTranscriptLocation().getStart().equals(stop)
                                                            && region.getTranscriptLocation().getStop().equals(start)) {
                                                        Location genomicLocation = new Location(a.getEnd(), a.getStart());
                                                        genomicLocation.setId(hearsayDAOBeanService.getLocationDAO().save(genomicLocation));
                                                        region.setRegionLocation(genomicLocation);
                                                        hearsayDAOBeanService.getRegionDAO().save(region);
                                                    }
                                                } else {
                                                    if (region.getTranscriptLocation().getStart().equals(start)
                                                            && region.getTranscriptLocation().getStop().equals(stop)) {
                                                        Location genomicLocation = new Location(a.getStart(), a.getEnd());
                                                        genomicLocation.setId(hearsayDAOBeanService.getLocationDAO().save(genomicLocation));
                                                        region.setRegionLocation(genomicLocation);
                                                        hearsayDAOBeanService.getRegionDAO().save(region);
                                                    }
                                                }
                                            }
                                        } catch (HearsayDAOException e) {
                                            e.printStackTrace();
                                        }

                                    });

                                }

                            } catch (Exception e) {
                                logger.error(e.getMessage(), e);
                                e.printStackTrace();
                            }

                        }

                    });
                }

            }
            es.shutdown();
            es.awaitTermination(4L, TimeUnit.HOURS);

        } catch (Exception e) {
            logger.error(e.getMessage(), e);
        }

    }

}
