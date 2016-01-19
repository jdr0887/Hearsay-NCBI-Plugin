package org.renci.hearsay.commands.ncbi;

import static org.renci.hearsay.commands.ncbi.Constants.IDENTIFIER_KEY_NUCCORE;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.collections4.CollectionUtils;
import org.apache.commons.lang3.Range;
import org.renci.hearsay.dao.HearsayDAOBeanService;
import org.renci.hearsay.dao.HearsayDAOException;
import org.renci.hearsay.dao.model.Alignment;
import org.renci.hearsay.dao.model.Location;
import org.renci.hearsay.dao.model.ReferenceSequence;
import org.renci.hearsay.dao.model.Region;
import org.renci.hearsay.dao.model.RegionType;
import org.renci.hearsay.dao.model.StrandType;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class AddAlignmentUTRsRunnable implements Runnable {

    private static final Logger logger = LoggerFactory.getLogger(AddAlignmentUTRsRunnable.class);

    private HearsayDAOBeanService hearsayDAOBeanService;

    public AddAlignmentUTRsRunnable(HearsayDAOBeanService hearsayDAOBeanService) {
        super();
        this.hearsayDAOBeanService = hearsayDAOBeanService;
    }

    @Override
    public void run() {
        logger.debug("ENTERING run()");

        try {

            List<ReferenceSequence> referenceSequences = hearsayDAOBeanService.getReferenceSequenceDAO()
                    .findByIdentifierSystem("includeAlignments", IDENTIFIER_KEY_NUCCORE);

            if (CollectionUtils.isEmpty(referenceSequences)) {
                logger.warn("no reference sequences found");
                return;
            }

            ExecutorService es = Executors.newFixedThreadPool(4);

            for (ReferenceSequence referenceSequence : referenceSequences) {

                es.submit(() -> {

                    logger.info(referenceSequence.toString());

                    try {

                        StrandType strandType = referenceSequence.getStrandType();

                        Set<Alignment> alignments = referenceSequence.getAlignments();

                        if (CollectionUtils.isEmpty(alignments)) {
                            logger.warn("no alignments found: {}", referenceSequence.toString());
                            return;
                        }

                        logger.info("alignments.size(): {}", alignments.size());

                        for (Alignment alignment : alignments) {

                            Location proteinLocation = alignment.getProteinLocation();
                            logger.info("Protein: {}", proteinLocation.toString());

                            List<Region> regionList = hearsayDAOBeanService.getRegionDAO().findByAlignmentId(alignment.getId());

                            if (CollectionUtils.isEmpty(regionList)) {
                                logger.warn("no regions found: {}", alignment.toString());
                                continue;
                            }

                            List<Region> utrRegionList = new ArrayList<Region>();

                            regionList.sort((a, b) -> {
                                if (a != null && a.getRegionLocation() != null && b != null && b.getRegionLocation() != null) {
                                    return a.getRegionLocation().getStart().compareTo(b.getRegionLocation().getStart());
                                }
                                return 0;
                            });

                            for (Region region : regionList) {
                                Location transcriptLocation = region.getTranscriptLocation();
                                if (transcriptLocation == null) {
                                    continue;
                                }
                                int transcriptStart = transcriptLocation.getStart();
                                int transcriptStop = transcriptLocation.getStop();
                                Range<Integer> transcriptRange = transcriptLocation.toRange();

                                Location regionLocation = region.getRegionLocation();
                                if (regionLocation == null) {
                                    continue;
                                }
                                int regionStart = regionLocation.getStart();
                                int regionStop = regionLocation.getStop();

                                if (strandType.equals(StrandType.MINUS)) {
                                    if (transcriptRange.contains(proteinLocation.getStart())
                                            && transcriptRange.contains(proteinLocation.getStop())) {

                                        transcriptLocation.setStop(proteinLocation.getStart() - 1);
                                        regionLocation.setStop(regionStart + transcriptLocation.diff());

                                        Region newRegion = createRegion(alignment, proteinLocation.getStart(), proteinLocation.getStop(),
                                                regionLocation.getStop() + 1,
                                                regionLocation.getStop() + 1 + (proteinLocation.getStop() - proteinLocation.getStart()));
                                        utrRegionList.add(newRegion);

                                        newRegion = createRegion(alignment, proteinLocation.getStop() + 1, transcriptStop,
                                                newRegion.getRegionLocation().getStop() + 1, regionStop);
                                        utrRegionList.add(newRegion);

                                    } else if (transcriptRange.contains(proteinLocation.getStart())) {

                                        transcriptLocation.setStop(proteinLocation.getStart() - 1);
                                        regionLocation.setStart(regionStop - transcriptLocation.diff());

                                        Region newRegion = createRegion(alignment, proteinLocation.getStart(), transcriptStop,
                                                regionLocation.getStart() - 1 - (transcriptStop - proteinLocation.getStart()),
                                                regionLocation.getStart() - 1);
                                        utrRegionList.add(newRegion);

                                    } else if (transcriptRange.contains(proteinLocation.getStop())) {

                                        transcriptLocation.setStart(proteinLocation.getStop() + 1);
                                        regionLocation.setStop(regionStart + transcriptLocation.diff());

                                        Region newRegion = createRegion(alignment, transcriptStart, proteinLocation.getStop(),
                                                regionLocation.getStop() + 1,
                                                regionLocation.getStop() + 1 + (proteinLocation.getStop() - transcriptStart));
                                        utrRegionList.add(newRegion);

                                    }

                                } else {
                                    if (transcriptRange.contains(proteinLocation.getStart())
                                            && transcriptRange.contains(proteinLocation.getStop())) {

                                        transcriptLocation.setStop(proteinLocation.getStart() - 1);
                                        regionLocation.setStop(regionStart + transcriptLocation.diff());

                                        Region newRegion = createRegion(alignment, proteinLocation.getStart(), proteinLocation.getStop(),
                                                regionLocation.getStop() + 1,
                                                regionLocation.getStop() + 1 + (proteinLocation.getStop() - proteinLocation.getStart()));
                                        utrRegionList.add(newRegion);

                                        newRegion = createRegion(alignment, proteinLocation.getStop() + 1, transcriptStop,
                                                newRegion.getRegionLocation().getStop() + 1, newRegion.getRegionLocation().getStop() + 1
                                                        + (transcriptStop - proteinLocation.getStop() - 1));
                                        utrRegionList.add(newRegion);

                                    } else if (transcriptRange.contains(proteinLocation.getStart())) {

                                        transcriptLocation.setStop(proteinLocation.getStart() - 1);
                                        regionLocation.setStop(regionStart + transcriptLocation.diff());

                                        Region newRegion = createRegion(alignment, proteinLocation.getStart(), transcriptStop,
                                                regionLocation.getStop() + 1,
                                                regionLocation.getStop() + 1 + (transcriptStop - proteinLocation.getStart()));
                                        utrRegionList.add(newRegion);

                                    } else if (transcriptRange.contains(proteinLocation.getStop())) {

                                        transcriptLocation.setStart(proteinLocation.getStop() + 1);
                                        regionLocation.setStart(regionStop - transcriptLocation.diff());

                                        Region newRegion = createRegion(alignment, transcriptStart, proteinLocation.getStop(),
                                                regionLocation.getStart() - 1 - (proteinLocation.getStop() - transcriptStart),
                                                regionLocation.getStart() - 1);
                                        utrRegionList.add(newRegion);

                                    }

                                }

                                hearsayDAOBeanService.getLocationDAO().save(transcriptLocation);
                                hearsayDAOBeanService.getLocationDAO().save(regionLocation);

                            }

                            regionList.addAll(utrRegionList);

                            for (Region region : regionList) {

                                if (strandType.equals(StrandType.PLUS)
                                        && region.getTranscriptLocation().getStop() < alignment.getProteinLocation().getStart()) {
                                    region.setRegionType(RegionType.UTR5);
                                }

                                if (strandType.equals(StrandType.PLUS)
                                        && region.getTranscriptLocation().getStop() > alignment.getProteinLocation().getStop()) {
                                    region.setRegionType(RegionType.UTR3);
                                }

                                if (strandType.equals(StrandType.MINUS)
                                        && region.getTranscriptLocation().getStop() < alignment.getProteinLocation().getStart()) {
                                    region.setRegionType(RegionType.UTR5);
                                }

                                if (strandType.equals(StrandType.MINUS)
                                        && region.getTranscriptLocation().getStop() > alignment.getProteinLocation().getStop()) {
                                    region.setRegionType(RegionType.UTR3);
                                }
                                hearsayDAOBeanService.getRegionDAO().save(region);

                            }

                            // adding intron regions

                            Region lastRegion = regionList.get(regionList.size() - 1);
                            Region previousRegion = null;
                            for (Region currentRegion : regionList) {
                                if (previousRegion == null) {
                                    previousRegion = currentRegion;
                                    continue;
                                }

                                if (previousRegion.getRegionLocation() == null || currentRegion.getRegionLocation() == null) {
                                    previousRegion = currentRegion;
                                    continue;
                                }

                                if (previousRegion.getRegionLocation().getStop().equals(currentRegion.getRegionLocation().getStart() - 1)) {
                                    previousRegion = currentRegion;
                                    continue;
                                }

                                if (currentRegion.equals(lastRegion)) {
                                    break;
                                }

                                Region region = new Region(RegionType.INTRON);
                                region.setAlignment(alignment);
                                Location regionLocation = new Location(previousRegion.getRegionLocation().getStop() + 1,
                                        currentRegion.getRegionLocation().getStart() - 1);
                                regionLocation.setId(hearsayDAOBeanService.getLocationDAO().save(regionLocation));
                                region.setRegionLocation(regionLocation);
                                hearsayDAOBeanService.getRegionDAO().save(region);

                                previousRegion = currentRegion;
                            }

                        }
                    } catch (Exception e) {
                        logger.error(e.getMessage(), e);
                        e.printStackTrace();
                    }

                });

            }
            es.shutdown();
            es.awaitTermination(2, TimeUnit.HOURS);
        } catch (Exception e) {
            logger.error(e.getMessage(), e);
            e.printStackTrace();
        }

    }

    private Region createRegion(Alignment alignment, Integer transcriptStart, Integer transcriptStop, Integer genomicStart,
            Integer genomicStop) throws HearsayDAOException {
        Region newRegion = new Region(RegionType.EXON);
        newRegion.setAlignment(alignment);

        Location newTranscriptLocation = new Location(transcriptStart, transcriptStop);
        newTranscriptLocation.setId(hearsayDAOBeanService.getLocationDAO().save(newTranscriptLocation));
        newRegion.setTranscriptLocation(newTranscriptLocation);

        Location newRegionLocation = new Location(genomicStart, genomicStop);
        newRegionLocation.setId(hearsayDAOBeanService.getLocationDAO().save(newRegionLocation));
        newRegion.setRegionLocation(newRegionLocation);

        newRegion.setId(hearsayDAOBeanService.getRegionDAO().save(newRegion));
        return newRegion;
    }

}
