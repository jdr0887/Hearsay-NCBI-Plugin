package org.renci.hearsay.commands.ncbi;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import org.apache.commons.lang3.math.NumberUtils;
import org.renci.gbff.model.Feature;
import org.renci.gbff.model.Sequence;
import org.renci.hearsay.commands.ncbi.util.FTPUtil;
import org.renci.hearsay.dao.HearsayDAOBean;
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

public class PersistAlignmentsFromUCSCRunnable implements Runnable {

    private static final Logger logger = LoggerFactory.getLogger(PersistAlignmentsFromUCSCRunnable.class);

    private static final LinkedList<String> refGeneLines = new LinkedList<String>();

    static {
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
    }

    private final Pattern featureLocationPattern = Pattern.compile("^(join|order)\\((.+)\\)$");

    private final HearsayDAOBean hearsayDAOBean;

    private final Sequence sequence;

    public PersistAlignmentsFromUCSCRunnable(final HearsayDAOBean hearsayDAOBean, final Sequence sequence) {
        super();
        this.hearsayDAOBean = hearsayDAOBean;
        this.sequence = sequence;
    }

    @Override
    public void run() {
        logger.debug("ENTERING run()");

        try {
            String refseqAccession = sequence.getAccession().contains(" ") ? sequence.getAccession().substring(0,
                    sequence.getAccession().indexOf(" ")) : sequence.getAccession();

            String refSeqVersionedAccession = sequence.getVersion().trim().contains(" ") ? sequence.getVersion()
                    .substring(0, sequence.getVersion().indexOf(" ")) : sequence.getVersion();

            List<Identifier> rnaNucleotideAccessionIdentifierList = hearsayDAOBean.getIdentifierDAO().findByExample(
                    new Identifier("www.ncbi.nlm.nih.gov/nuccore", refSeqVersionedAccession));

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

            List<Identifier> proteinAccessionIdentifierList = hearsayDAOBean.getIdentifierDAO().findByExample(
                    new Identifier("www.ncbi.nlm.nih.gov/protein", proteinAccession));

            List<Long> identifierIdList = new ArrayList<Long>();
            for (Identifier identifier : rnaNucleotideAccessionIdentifierList) {
                identifierIdList.add(identifier.getId());
            }
            for (Identifier identifier : proteinAccessionIdentifierList) {
                identifierIdList.add(identifier.getId());
            }

            List<ReferenceSequence> potentialRefSeqs = hearsayDAOBean.getReferenceSequenceDAO().findByIdentifiers(
                    identifierIdList);

            if (potentialRefSeqs == null || (potentialRefSeqs != null && potentialRefSeqs.isEmpty())) {
                logger.warn("Could not find ReferenceSequence: refSeqVersionedAccession = {}, proteinAccession = {}",
                        refSeqVersionedAccession, proteinAccession);
                return;
            }

            logger.info("Using ReferenceSequence: refSeqVersionedAccession = {}, proteinAccession = {}",
                    refSeqVersionedAccession, proteinAccession);

            final StrandType strandType = potentialRefSeqs.get(0).getStrandType();

            // add alignments
            Alignment alignment = new Alignment();
            alignment.setId(hearsayDAOBean.getAlignmentDAO().save(alignment));
            alignment.setReferenceSequences(potentialRefSeqs);

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
            proteinLocation.setId(hearsayDAOBean.getLocationDAO().save(proteinLocation));
            alignment.setProteinLocation(proteinLocation);

            hearsayDAOBean.getAlignmentDAO().save(alignment);

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
                    transcriptLocation.setId(hearsayDAOBean.getLocationDAO().save(transcriptLocation));
                    logger.debug(transcriptLocation.toString());
                    region.setTranscriptLocation(transcriptLocation);
                }
                regions.add(region);
            }

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
                        Location genomicLocation = new Location(Integer.valueOf(exonStartsSplit[i]),
                                Integer.valueOf(exonEndsSplit[i]));
                        logger.debug(genomicLocation.toString());
                        genomicLocationList.add(genomicLocation);
                    }
                    scanner.close();
                }
            }

            if (regions.size() != genomicLocationList.size()) {
                logger.warn("regions.size() = {}, genomicLocationList.size() = {}", regions.size(),
                        genomicLocationList.size());
                return;
            }

            Collections.sort(genomicLocationList, new Comparator<Location>() {
                @Override
                public int compare(Location l1, Location l2) {
                    return l1.getStart().compareTo(l2.getStart());
                }
            });

            // initial exon regions
            for (int i = 0; i < regions.size(); i++) {
                Region region = regions.get(i);
                region.setAlignment(alignment);
                Location regionLocation = genomicLocationList.get(i);
                regionLocation.setId(hearsayDAOBean.getLocationDAO().save(regionLocation));
                region.setRegionLocation(regionLocation);
                region.setId(hearsayDAOBean.getRegionDAO().save(region));
                alignment.getRegions().add(region);
            }

            hearsayDAOBean.getAlignmentDAO().save(alignment);

        } catch (Exception e) {
            logger.error(e.getMessage(), e);
            e.printStackTrace();
        }

    }

}
