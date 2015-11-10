package org.renci.hearsay.commands.ncbi;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.collections4.CollectionUtils;
import org.apache.commons.lang3.math.NumberUtils;
import org.renci.gbff.model.Feature;
import org.renci.gbff.model.Sequence;
import org.renci.hearsay.dao.HearsayDAOBean;
import org.renci.hearsay.dao.model.Alignment;
import org.renci.hearsay.dao.model.Identifier;
import org.renci.hearsay.dao.model.Location;
import org.renci.hearsay.dao.model.ReferenceSequence;
import org.renci.hearsay.dao.model.Region;
import org.renci.hearsay.dao.model.RegionType;
import org.renci.hearsay.dao.model.StrandType;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class PersistAlignmentsFromNCBIRunnable implements Runnable {

    private static final Logger logger = LoggerFactory.getLogger(PersistAlignmentsFromNCBIRunnable.class);

    private static final Pattern featureLocationPattern = Pattern.compile("^(join|order)\\((.+)\\)$");

    private final HearsayDAOBean hearsayDAOBean;

    private final Sequence sequence;

    private final List<ReferenceSequence> referenceSequences;

    private final Feature firstCDSFeature;

    private final LinkedList<String> alignmentsLines;

    public PersistAlignmentsFromNCBIRunnable(final HearsayDAOBean hearsayDAOBean, final Sequence sequence,
            final List<ReferenceSequence> referenceSequences, final Feature firstCDSFeature,
            final LinkedList<String> alignmentsLines) {
        super();
        this.hearsayDAOBean = hearsayDAOBean;
        this.sequence = sequence;
        this.referenceSequences = referenceSequences;
        this.firstCDSFeature = firstCDSFeature;
        this.alignmentsLines = alignmentsLines;
    }

    @Override
    public void run() {
        logger.debug("ENTERING run()");

        try {
            String refSeqVersionedAccession = sequence.getVersion().trim().contains(" ") ? sequence.getVersion()
                    .substring(0, sequence.getVersion().indexOf(" ")) : sequence.getVersion();

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
            alignment.setId(hearsayDAOBean.getAlignmentDAO().save(alignment));
            alignment.setReferenceSequences(referenceSequences);

            Location proteinLocation = new Location(firstCDSFeatureLocationStart, firstCDSFeatureLocationStop);
            proteinLocation.setId(hearsayDAOBean.getLocationDAO().save(proteinLocation));
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
                    transcriptLocation.setId(hearsayDAOBean.getLocationDAO().save(transcriptLocation));
                    region.setTranscriptLocation(transcriptLocation);
                }

                region.setId(hearsayDAOBean.getRegionDAO().save(region));
                alignment.getRegions().add(region);
            }

            hearsayDAOBean.getAlignmentDAO().save(alignment);

            List<Region> regions = alignment.getRegions();
            if (CollectionUtils.isNotEmpty(regions)) {
                Collections.sort(regions, new Comparator<Region>() {
                    @Override
                    public int compare(Region r1, Region r2) {
                        if (strandType.equals(StrandType.MINUS)) {
                            return r2.getTranscriptLocation().getStart()
                                    .compareTo(r1.getTranscriptLocation().getStart());
                        } else {
                            return r1.getTranscriptLocation().getStart()
                                    .compareTo(r2.getTranscriptLocation().getStart());
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
                                    genomicLocation.setId(hearsayDAOBean.getLocationDAO().save(genomicLocation));
                                    region.setRegionLocation(genomicLocation);
                                    hearsayDAOBean.getRegionDAO().save(region);
                                }
                            } else {
                                if (region.getTranscriptLocation().getStart().equals(start)
                                        && region.getTranscriptLocation().getStop().equals(stop)) {
                                    Location genomicLocation = new Location(Integer.valueOf(genomicStart),
                                            Integer.valueOf(genomicStop));
                                    genomicLocation.setId(hearsayDAOBean.getLocationDAO().save(genomicLocation));
                                    region.setRegionLocation(genomicLocation);
                                    hearsayDAOBean.getRegionDAO().save(region);
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
