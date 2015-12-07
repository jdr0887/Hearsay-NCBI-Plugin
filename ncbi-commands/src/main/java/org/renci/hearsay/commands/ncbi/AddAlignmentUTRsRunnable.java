package org.renci.hearsay.commands.ncbi;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
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

            ExecutorService es = Executors.newFixedThreadPool(4);
            List<Alignment> alignmentList = hearsayDAOBeanService.getAlignmentDAO().findAll();
            if (CollectionUtils.isNotEmpty(alignmentList)) {
                logger.debug("alignmentList.size(): {}", alignmentList.size());
                for (Alignment alignment : alignmentList) {
                    if (CollectionUtils.isNotEmpty(alignment.getRegions())) {
                        es.submit(new AddUTRRunnable(alignment));
                    }
                }
            }
            es.shutdown();
            es.awaitTermination(2, TimeUnit.HOURS);
        } catch (Exception e) {
            logger.error(e.getMessage(), e);
            e.printStackTrace();
        }

    }

    class AddUTRRunnable implements Runnable {

        private Alignment alignment;

        public AddUTRRunnable(Alignment alignment) {
            super();
            this.alignment = alignment;
        }

        @Override
        public void run() {
            logger.debug("ENTERING run()");

            try {

                List<Region> regionList = alignment.getRegions();
                Set<ReferenceSequence> referenceSequenceSet = alignment.getReferenceSequences();

                if (CollectionUtils.isEmpty(referenceSequenceSet)) {
                    logger.error("Could not find ReferenceSequence for Alignment: {}", alignment.toString());
                    return;
                }

                StrandType strandType = referenceSequenceSet.iterator().next().getStrandType();
                List<Region> utrRegionList = new ArrayList<Region>();

                final Location proteinLocation = alignment.getProteinLocation();
                logger.info("Protein: {}", proteinLocation.toString());

                // regionList.sort((a, b) ->
                // a.getRegionLocation().getStart().compareTo(b.getRegionLocation().getStart()));

                Collections.sort(regionList, new Comparator<Region>() {
                    @Override
                    public int compare(Region o1, Region o2) {
                        Location first = o1.getRegionLocation();
                        Location second = o2.getRegionLocation();
                        if (first != null && second != null) {
                            return Integer.compare(first.getStart(), second.getStart());
                        }
                        return 0;
                    }
                });

                for (Region region : regionList) {
                    Location regionLocation = region.getRegionLocation();
                    if (regionLocation == null) {
                        continue;
                    }
                    int regionStart = regionLocation.getStart();
                    int regionStop = regionLocation.getStop();

                    Location transcriptLocation = region.getTranscriptLocation();
                    if (transcriptLocation == null) {
                        continue;
                    }
                    int transcriptStart = transcriptLocation.getStart();
                    int transcriptStop = transcriptLocation.getStop();
                    Range<Integer> transcriptRange = transcriptLocation.toRange();

                    if (strandType.equals(StrandType.MINUS)) {

                        if (transcriptRange.contains(proteinLocation.getStart())
                                && transcriptRange.contains(proteinLocation.getStop())) {

                            transcriptLocation.setStop(proteinLocation.getStart() - 1);
                            regionLocation.setStop(regionStop - transcriptLocation.diff());

                            Region newRegion = createRegion(alignment, proteinLocation.getStart(),
                                    proteinLocation.getStop(), regionLocation.getStop() + 1, regionLocation.getStop()
                                            + 1 + (proteinLocation.getStop() - proteinLocation.getStart()));
                            utrRegionList.add(newRegion);

                            newRegion = createRegion(alignment, proteinLocation.getStop() + 1, transcriptStop,
                                    newRegion.getRegionLocation().getStop() + 1, newRegion.getRegionLocation().getStop()
                                            + 1 + (transcriptStop - proteinLocation.getStop() - 1));
                            utrRegionList.add(newRegion);

                        } else if (transcriptRange.contains(proteinLocation.getStart())) {

                            transcriptLocation.setStart(proteinLocation.getStart() - 1);
                            regionLocation.setStart(regionStop - transcriptLocation.diff());

                            Region newRegion = createRegion(alignment, transcriptStart, proteinLocation.getStart(),
                                    regionLocation.getStart() - 1 - (transcriptStart - proteinLocation.getStart()),
                                    regionLocation.getStart() - 1);
                            utrRegionList.add(newRegion);

                        } else if (transcriptRange.contains(proteinLocation.getStop())) {

                            transcriptLocation.setStart(proteinLocation.getStop());
                            regionLocation.setStart(regionLocation.getStop() - transcriptLocation.diff());

                            Region newRegion = createRegion(alignment, transcriptStart, proteinLocation.getStop() + 1,
                                    regionLocation.getStart() - (transcriptStart - proteinLocation.getStop()),
                                    regionLocation.getStart() - 1);
                            utrRegionList.add(newRegion);

                        }

                    }

                    if (strandType.equals(StrandType.PLUS)) {

                        if (transcriptRange.contains(proteinLocation.getStart())
                                && transcriptRange.contains(proteinLocation.getStop())) {

                            transcriptLocation.setStop(proteinLocation.getStart() - 1);
                            regionLocation.setStop(regionStart + transcriptLocation.diff());

                            Region newRegion = createRegion(alignment, proteinLocation.getStart(),
                                    proteinLocation.getStop(), regionLocation.getStop() + 1, regionLocation.getStop()
                                            + 1 + (proteinLocation.getStop() - proteinLocation.getStart()));
                            utrRegionList.add(newRegion);

                            newRegion = createRegion(alignment, proteinLocation.getStop() + 1, transcriptStop,
                                    newRegion.getRegionLocation().getStop() + 1, newRegion.getRegionLocation().getStop()
                                            + 1 + (transcriptStop - proteinLocation.getStop() - 1));
                            utrRegionList.add(newRegion);

                        } else if (transcriptRange.contains(proteinLocation.getStart())) {

                            transcriptLocation.setStop(proteinLocation.getStart() - 1);
                            regionLocation.setStop(regionStart + transcriptLocation.diff());

                            Region newRegion = createRegion(alignment, proteinLocation.getStart(), transcriptStop,
                                    regionLocation.getStop() + 1,
                                    regionLocation.getStop() + 1 + (transcriptStop - proteinLocation.getStart()));
                            utrRegionList.add(newRegion);

                        } else if (transcriptRange.contains(proteinLocation.getStop())) {

                            transcriptLocation.setStop(proteinLocation.getStop() - 1);
                            regionLocation.setStop(regionStart + transcriptLocation.diff());

                            Region newRegion = createRegion(alignment, proteinLocation.getStop(), transcriptStop,
                                    regionStop - (transcriptStop - proteinLocation.getStop() + 1), regionStop);
                            utrRegionList.add(newRegion);

                        }

                    }
                }

                regionList.addAll(utrRegionList);

                Collections.sort(regionList, new Comparator<Region>() {
                    @Override
                    public int compare(Region o1, Region o2) {
                        Location first = o1.getRegionLocation();
                        Location second = o2.getRegionLocation();
                        if (first != null && second != null) {
                            return Integer.compare(first.getStart(), second.getStart());
                        }
                        return 0;
                    }
                });

                for (Region region : regionList) {
                    Location transcriptLocation = region.getTranscriptLocation();
                    if (transcriptLocation == null) {
                        continue;
                    }

                    if (strandType.equals(StrandType.PLUS)
                            && transcriptLocation.getStop() < proteinLocation.getStart()) {
                        region.setRegionType(RegionType.UTR5);
                    }

                    if (strandType.equals(StrandType.PLUS)
                            && transcriptLocation.getStop() > proteinLocation.getStop()) {
                        region.setRegionType(RegionType.UTR3);
                    }

                    if (strandType.equals(StrandType.MINUS)
                            && transcriptLocation.getStop() < proteinLocation.getStart()) {
                        region.setRegionType(RegionType.UTR5);
                    }

                    if (strandType.equals(StrandType.MINUS)
                            && transcriptLocation.getStop() > proteinLocation.getStop()) {
                        region.setRegionType(RegionType.UTR3);
                    }

                    hearsayDAOBeanService.getRegionDAO().save(region);
                }

                // for (Region region : regionList) {
                // Location regionLocation = region.getRegionLocation();
                // Location transcriptLocation = region.getTranscriptLocation();
                // if (transcriptLocation == null) {
                // logger.info("{}, {}", region.toString(),
                // regionLocation.toString());
                // continue;
                // }
                // logger.info("{}, {}, {}", region.toString(),
                // regionLocation.toString(),
                // transcriptLocation.toString());
                // }

                List<Region> intronRegionList = new ArrayList<Region>();
                Region lastRegion = regionList.get(regionList.size() - 1);

                // adding intron regions
                Region previousRegion = null;
                for (Region currentRegion : regionList) {
                    if (previousRegion == null) {
                        previousRegion = currentRegion;
                        continue;
                    }

                    Location previousRegionLocation = previousRegion.getRegionLocation();

                    if (previousRegionLocation == null) {
                        continue;
                    }

                    Location currentRegionLocation = currentRegion.getRegionLocation();

                    if (currentRegionLocation == null) {
                        continue;
                    }

                    if (previousRegionLocation.getStop().equals(currentRegionLocation.getStart() - 1)) {
                        previousRegion = currentRegion;
                        continue;
                    }

                    if (currentRegion.equals(lastRegion)) {
                        break;
                    }

                    Region region = new Region(RegionType.INTRON);
                    region.setId(hearsayDAOBeanService.getRegionDAO().save(region));
                    region.setAlignment(alignment);

                    Location regionLocation = new Location(previousRegionLocation.getStop() + 1,
                            currentRegionLocation.getStart() - 1);
                    regionLocation.setId(hearsayDAOBeanService.getLocationDAO().save(regionLocation));
                    region.setRegionLocation(regionLocation);
                    hearsayDAOBeanService.getRegionDAO().save(region);
                    intronRegionList.add(region);

                    previousRegion = currentRegion;
                }

                regionList.addAll(intronRegionList);

                Collections.sort(regionList, new Comparator<Region>() {
                    @Override
                    public int compare(Region o1, Region o2) {
                        Location first = o1.getRegionLocation();
                        Location second = o2.getRegionLocation();
                        if (first != null && second != null) {
                            return Integer.compare(first.getStart(), second.getStart());
                        }
                        return 0;
                    }
                });

                for (Region region : regionList) {
                    Location regionLocation = region.getRegionLocation();
                    if (regionLocation == null) {
                        continue;
                    }
                    Location transcriptLocation = region.getTranscriptLocation();
                    if (transcriptLocation == null) {
                        logger.debug("{}, {}", region.toString(), regionLocation.toString());
                        continue;
                    }
                    logger.debug("{}, {}, {}", region.toString(), regionLocation.toString(),
                            transcriptLocation.toString());
                }

                alignment.setRegions(regionList);

                hearsayDAOBeanService.getAlignmentDAO().save(alignment);

            } catch (Exception e) {
                logger.error("Error", e);
                e.printStackTrace();
            }

        }

        private Region createRegion(Alignment alignment, Integer transcriptStart, Integer transcriptStop,
                Integer genomicStart, Integer genomicStop) throws HearsayDAOException {
            Region newRegion = new Region(RegionType.EXON);
            newRegion.setId(hearsayDAOBeanService.getRegionDAO().save(newRegion));
            newRegion.setAlignment(alignment);

            Location newTranscriptLocation = new Location(transcriptStart, transcriptStop);
            newTranscriptLocation.setId(hearsayDAOBeanService.getLocationDAO().save(newTranscriptLocation));
            newRegion.setTranscriptLocation(newTranscriptLocation);

            Location newRegionLocation = new Location(genomicStart, genomicStop);
            newRegionLocation.setId(hearsayDAOBeanService.getLocationDAO().save(newRegionLocation));
            newRegion.setRegionLocation(newRegionLocation);
            hearsayDAOBeanService.getRegionDAO().save(newRegion);
            return newRegion;
        }
    }

}
