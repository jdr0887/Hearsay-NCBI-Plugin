package org.renci.hearsay.commands.ncbi;

import static org.renci.hearsay.commands.ncbi.Constants.IDENTIFIER_KEY_NUCCORE;
import static org.renci.hearsay.commands.ncbi.Constants.IDENTIFIER_KEY_PROTEIN;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.collections4.CollectionUtils;
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
import org.renci.hearsay.dao.model.Alignment;
import org.renci.hearsay.dao.model.Identifier;
import org.renci.hearsay.dao.model.Location;
import org.renci.hearsay.dao.model.ReferenceSequence;
import org.renci.hearsay.dao.model.Region;
import org.renci.hearsay.dao.model.RegionType;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class PullAlignmentsRunnable implements Runnable {

    private static final Logger logger = LoggerFactory.getLogger(PullAlignmentsRunnable.class);

    private static final Pattern featureLocationPattern = Pattern.compile("^(join|order)\\((.+)\\)$");

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

            fileList.forEach(a -> logger.info(a.getAbsolutePath()));

            for (File f : fileList) {

                logger.info("parsing GenBankFlatFile: {}", f.getAbsolutePath());

                List<Sequence> sequenceList = gbffMgr.deserialize(gbffFilter, f);

                if (CollectionUtils.isEmpty(sequenceList)) {
                    logger.warn("no sequences found");
                    continue;
                }

                logger.info("sequenceList.size(): {}", sequenceList.size());

                ExecutorService es = Executors.newFixedThreadPool(4);

                for (Sequence sequence : sequenceList) {

                    es.submit(() -> {

                        try {

                            logger.info(sequence.toString());

                            if (CollectionUtils.isEmpty(sequence.getFeatures())) {
                                logger.warn("sequence.getFeatures() is empty");
                                return;
                            }

                            List<Identifier> identifierList = new ArrayList<Identifier>();

                            // rna nucleotide accession
                            String refSeqVersionedAccession = sequence.getVersion().trim().contains(" ")
                                    ? sequence.getVersion().substring(0, sequence.getVersion().indexOf(" ")) : sequence.getVersion();

                            List<Identifier> rnaNucleotideAccessionIdentifierList = hearsayDAOBeanService.getIdentifierDAO()
                                    .findByExample(new Identifier(IDENTIFIER_KEY_NUCCORE, refSeqVersionedAccession));
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
                            proteinAccession = firstCDSFeature.getQualifiers().get("protein_id").replace("\"", "");

                            int exonCount = 0;
                            for (Feature feature : sequence.getFeatures()) {
                                if ("exon".equals(feature.getType())) {
                                    exonCount++;
                                }
                            }

                            if (exonCount == 0) {
                                logger.warn("no exons found: {}", sequence.toString());
                                return;
                            }

                            logger.info("number of exons found: {}", exonCount);

                            List<Identifier> proteinAccessionIdentifierList = hearsayDAOBeanService.getIdentifierDAO()
                                    .findByExample(new Identifier(IDENTIFIER_KEY_PROTEIN, proteinAccession));
                            if (CollectionUtils.isNotEmpty(proteinAccessionIdentifierList)) {
                                identifierList.add(proteinAccessionIdentifierList.get(0));
                            }

                            identifierList.forEach(a -> logger.info(a.toString()));

                            List<ReferenceSequence> potentialRefSeqs = hearsayDAOBeanService.getReferenceSequenceDAO()
                                    .findByIdentifiers(identifierList);

                            if (CollectionUtils.isEmpty(potentialRefSeqs)) {
                                logger.warn("Could not find ReferenceSequence");
                                return;
                            }

                            List<GFF3Record> gff3Records = gff3Mgr
                                    .deserialize(new AttributeValueFilter("Target", refSeqVersionedAccession));
                            logger.info("gff3Records.size(): {}", gff3Records.size());
                            if (CollectionUtils.isEmpty(gff3Records)) {
                                logger.warn("gff3Records is empty");
                                return;
                            }

                            ReferenceSequence referenceSequence = potentialRefSeqs.get(0);
                            logger.info(referenceSequence.toString());

                            // add protein info to alignment
                            Location proteinLocation = null;
                            String firstCDSFeatureLocation = firstCDSFeature.getLocation();
                            logger.debug("firstCDSFeatureLocation: {}", firstCDSFeatureLocation);
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
                                logger.info("proteinLocation: {}", proteinLocation.toString());
                                alignment.setProteinLocation(proteinLocation);
                            }
                            alignment.getReferenceSequences().add(referenceSequence);
                            alignment.setId(hearsayDAOBeanService.getAlignmentDAO().save(alignment));

                            // add exons to alignment
                            for (Feature feature : sequence.getFeatures()) {

                                if (!"exon".equals(feature.getType())) {
                                    continue;
                                }

                                try {

                                    String range = feature.getLocation();
                                    String[] split = range.split("\\.\\.");

                                    Location transcriptLocation = null;
                                    if (NumberUtils.isNumber(split[0]) && NumberUtils.isNumber(split[1])) {
                                        transcriptLocation = new Location(Integer.valueOf(split[0]), Integer.valueOf(split[1]));
                                        transcriptLocation.setId(hearsayDAOBeanService.getLocationDAO().save(transcriptLocation));
                                        logger.debug("transcriptLocation: {}", transcriptLocation.toString());
                                    }

                                    if (transcriptLocation == null) {
                                        logger.warn("exon with null transcript: {}", sequence.toString());
                                        continue;
                                    }

                                    Location genomicLocation = null;
                                    for (GFF3Record record : gff3Records) {
                                        String targetValue = record.getAttributes().get("Target");
                                        String[] targetSplit = targetValue.split(" ");
                                        Integer start = Integer.valueOf(targetSplit[1]);
                                        Integer stop = Integer.valueOf(targetSplit[2]);
                                        if (transcriptLocation.getStart().equals(start) && transcriptLocation.getStop().equals(stop)) {
                                            genomicLocation = new Location(record.getStart(), record.getEnd());
                                            genomicLocation.setId(hearsayDAOBeanService.getLocationDAO().save(genomicLocation));
                                            logger.debug("genomicLocation: {}", genomicLocation.toString());
                                            break;
                                        }
                                    }
                                    Region region = new Region(RegionType.EXON);
                                    region.setTranscriptLocation(transcriptLocation);
                                    region.setRegionLocation(genomicLocation);
                                    region.setAlignment(alignment);
                                    region.setId(hearsayDAOBeanService.getRegionDAO().save(region));
                                } catch (Exception e) {
                                    logger.error(e.getMessage(), e);
                                    e.printStackTrace();
                                }

                            }

                            
                        } catch (Exception e) {
                            logger.error(e.getMessage(), e);
                        }
                    });
                }
                es.shutdown();
                es.awaitTermination(1L, TimeUnit.HOURS);
            }
        } catch (Exception e) {
            logger.error(e.getMessage(), e);
        }

    }

}
