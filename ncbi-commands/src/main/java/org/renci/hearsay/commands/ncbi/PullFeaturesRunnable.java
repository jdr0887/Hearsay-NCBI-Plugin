package org.renci.hearsay.commands.ncbi;

import static org.renci.hearsay.commands.ncbi.Constants.IDENTIFIER_KEY_NUCCORE;
import static org.renci.hearsay.commands.ncbi.Constants.IDENTIFIER_KEY_PROTEIN;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;
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
import org.renci.hearsay.commands.ncbi.util.FTPUtil;
import org.renci.hearsay.dao.HearsayDAOBeanService;
import org.renci.hearsay.dao.model.Identifier;
import org.renci.hearsay.dao.model.Location;
import org.renci.hearsay.dao.model.ReferenceSequence;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class PullFeaturesRunnable implements Runnable {

    private static final Logger logger = LoggerFactory.getLogger(PullFeaturesRunnable.class);

    private static final Pattern featureLocationPattern = Pattern.compile("^(join|order)\\((.+)\\)$");

    private static final List<String> inclusionPatterns = Arrays.asList(new String[] { "misc_feature", "polyA_signal", "polyA_site",
            "transit_peptide", "mat_peptide", "sig_peptide", "unsure", "stem_loop", "protein_bind", "repeat_region", "prim_transcript",
            "proprotein", "LTR", "TATA_signal", "primer_bind", "terminator", "misc_difference", "misc_binding", "RBS", "misc_signal",
            "J_segment", "C_region", "conflict", "promoter", "ncRNA", "modified_base" });

    private HearsayDAOBeanService hearsayDAOBeanService;

    public PullFeaturesRunnable(HearsayDAOBeanService hearsayDAOBeanService) {
        super();
        this.hearsayDAOBeanService = hearsayDAOBeanService;
    }

    @Override
    public void run() {
        logger.debug("ENTERING run()");

        try {

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
                    return;
                }

                logger.info("sequenceList.size(): {}", sequenceList.size());

                for (Sequence sequence : sequenceList) {
                    logger.info(sequence.toString());

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
                        continue;
                    }

                    ReferenceSequence referenceSequence = potentialRefSeqs.get(0);
                    logger.info(referenceSequence.toString());

                    // add features
                    for (Feature feature : sequence.getFeatures()) {
                        if (!inclusionPatterns.contains(feature.getType())) {
                            continue;
                        }
                        org.renci.hearsay.dao.model.Feature hearsayFeature = new org.renci.hearsay.dao.model.Feature(feature.getType());
                        String note = feature.getQualifiers().get("note");
                        if (StringUtils.isNotEmpty(note)) {
                            hearsayFeature.setNote(note);
                        }
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

                        hearsayFeature.setId(hearsayDAOBeanService.getFeatureDAO().save(hearsayFeature));
                        logger.info(hearsayFeature.toString());

                    }

                }
            }
        } catch (Exception e) {
            logger.error(e.getMessage(), e);
        }
        logger.debug("LEAVING run()");
    }

}
