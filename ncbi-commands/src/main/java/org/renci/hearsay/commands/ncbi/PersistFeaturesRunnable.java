package org.renci.hearsay.commands.ncbi;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang3.math.NumberUtils;
import org.renci.gbff.model.Feature;
import org.renci.gbff.model.Sequence;
import org.renci.hearsay.dao.HearsayDAOBean;
import org.renci.hearsay.dao.model.Identifier;
import org.renci.hearsay.dao.model.Location;
import org.renci.hearsay.dao.model.ReferenceSequence;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class PersistFeaturesRunnable implements Runnable {

    private final Logger logger = LoggerFactory.getLogger(PersistFeaturesRunnable.class);

    private final Pattern featureLocationPattern = Pattern.compile("^(join|order)\\((.+)\\)$");

    private final List<String> inclusionPatterns = Arrays.asList(new String[] { "misc_feature", "polyA_signal",
            "polyA_site", "transit_peptide", "mat_peptide", "sig_peptide", "unsure", "stem_loop", "protein_bind",
            "repeat_region", "prim_transcript", "proprotein", "LTR", "TATA_signal", "primer_bind", "terminator",
            "misc_difference", "misc_binding", "RBS", "misc_signal", "J_segment", "C_region", "conflict", "promoter",
            "ncRNA", "modified_base" });

    private final HearsayDAOBean hearsayDAOBean;

    private final Sequence sequence;

    public PersistFeaturesRunnable(final HearsayDAOBean hearsayDAOBean, final Sequence sequence) {
        super();
        this.hearsayDAOBean = hearsayDAOBean;
        this.sequence = sequence;
    }

    @Override
    public void run() {
        logger.debug("ENTERING run()");

        try {

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

            logger.info(potentialRefSeqs.get(0).toString());

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

                hearsayFeature.setId(hearsayDAOBean.getFeatureDAO().save(hearsayFeature));
                hearsayFeature.setReferenceSequences(potentialRefSeqs);

                if (NumberUtils.isNumber(location)) {
                    Location l = new Location(Integer.valueOf(location), Integer.valueOf(location));
                    l.setId(hearsayDAOBean.getLocationDAO().save(l));
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
                                l.setId(hearsayDAOBean.getLocationDAO().save(l));
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
                        l.setId(hearsayDAOBean.getLocationDAO().save(l));
                        hearsayFeature.getLocations().add(l);
                    }
                }

                hearsayDAOBean.getFeatureDAO().save(hearsayFeature);
                logger.info(hearsayFeature.toString());

            }
        } catch (Exception e) {
            logger.error(e.getMessage(), e);
            e.printStackTrace();
        }

    }

}
