package org.renci.hearsay.commands.ncbi;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.lang.math.IntRange;
import org.renci.hearsay.dao.HearsayDAOBean;
import org.renci.hearsay.dao.model.Alignment;
import org.renci.hearsay.dao.model.Gene;
import org.renci.hearsay.dao.model.Location;
import org.renci.hearsay.dao.model.ReferenceSequence;
import org.renci.hearsay.dao.model.Region;
import org.renci.hearsay.dao.model.RegionType;
import org.renci.hearsay.dao.model.StrandType;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class AddAlignmentUTRsRunnable implements Runnable {

    private final Logger logger = LoggerFactory.getLogger(AddAlignmentUTRsRunnable.class);

    private HearsayDAOBean hearsayDAOBean;

    public AddAlignmentUTRsRunnable() {
        super();
    }

    @Override
    public void run() {
        logger.debug("ENTERING run()");

        try {

            List<Gene> geneList = hearsayDAOBean.getGeneDAO().findAll();
            if (CollectionUtils.isNotEmpty(geneList)) {
                for (Gene gene : geneList) {

                    logger.info(gene.toString());

                    List<ReferenceSequence> referenceSequenceList = hearsayDAOBean.getReferenceSequenceDAO()
                            .findByGeneId(gene.getId());

                    if (CollectionUtils.isNotEmpty(referenceSequenceList)) {

                        for (ReferenceSequence referenceSequence : referenceSequenceList) {

                            StrandType strandType = referenceSequence.getStrandType();

                            logger.info(referenceSequence.toString());

                            List<Alignment> alignmentList = hearsayDAOBean.getAlignmentDAO().findByReferenceSequenceId(
                                    referenceSequence.getId());

                            if (CollectionUtils.isNotEmpty(alignmentList)) {

                                logger.debug("alignmentList.size(): {}", alignmentList.size());

                                for (Alignment alignment : alignmentList) {
                                    Location proteinLocation = alignment.getProteinLocation();
                                    logger.info("Protein: {}", proteinLocation.toString());

                                    List<Region> regionList = alignment.getRegions();

                                    List<Region> utrRegionList = new ArrayList<Region>();

                                    if (CollectionUtils.isNotEmpty(regionList)) {

                                        Collections.sort(regionList, new Comparator<Region>() {
                                            @Override
                                            public int compare(Region o1, Region o2) {
                                                return Integer.compare(o1.getRegionLocation().getStart(), o2
                                                        .getRegionLocation().getStart());
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
                                            IntRange transcriptRange = transcriptLocation.toRange();

                                            if (strandType.equals(StrandType.MINUS)) {

                                                if (transcriptRange.containsInteger(proteinLocation.getStart())) {

                                                    transcriptLocation.setStart(proteinLocation.getStart() - 1);
                                                    hearsayDAOBean.getLocationDAO().save(transcriptLocation);

                                                    regionLocation.setStart(regionLocation.getStop()
                                                            - transcriptLocation.diff());
                                                    hearsayDAOBean.getLocationDAO().save(regionLocation);

                                                    Region newRegion = new Region(RegionType.EXON);
                                                    newRegion.setId(hearsayDAOBean.getRegionDAO().save(newRegion));
                                                    newRegion.setAlignment(alignment);

                                                    Location newTranscriptLocation = new Location(transcriptStart,
                                                            proteinLocation.getStart());
                                                    newTranscriptLocation.setId(hearsayDAOBean.getLocationDAO().save(
                                                            newTranscriptLocation));
                                                    newRegion.setTranscriptLocation(newTranscriptLocation);

                                                    Location newRegionLocation = new Location(regionLocation.getStart()
                                                            - 1 - newTranscriptLocation.diff(),
                                                            regionLocation.getStart() - 1);
                                                    newRegionLocation.setId(hearsayDAOBean.getLocationDAO().save(
                                                            newRegionLocation));
                                                    newRegion.setRegionLocation(newRegionLocation);
                                                    hearsayDAOBean.getRegionDAO().save(newRegion);

                                                    utrRegionList.add(newRegion);
                                                }

                                                if (transcriptRange.containsInteger(proteinLocation.getStop())) {

                                                    transcriptLocation.setStart(proteinLocation.getStop());
                                                    hearsayDAOBean.getLocationDAO().save(transcriptLocation);

                                                    regionLocation.setStart(regionLocation.getStop()
                                                            - transcriptLocation.diff());
                                                    hearsayDAOBean.getLocationDAO().save(regionLocation);

                                                    Region newRegion = new Region(RegionType.EXON);
                                                    newRegion.setId(hearsayDAOBean.getRegionDAO().save(newRegion));
                                                    newRegion.setAlignment(alignment);

                                                    Location newTranscriptLocation = new Location(transcriptStart,
                                                            alignment.getProteinLocation().getStop() + 1);
                                                    newTranscriptLocation.setId(hearsayDAOBean.getLocationDAO().save(
                                                            newTranscriptLocation));
                                                    newRegion.setTranscriptLocation(newTranscriptLocation);

                                                    // is this off by one?
                                                    Location newRegionLocation = new Location(regionLocation.getStart()
                                                            - newTranscriptLocation.diff() - 1,
                                                            regionLocation.getStart() - 1);
                                                    newRegionLocation.setId(hearsayDAOBean.getLocationDAO().save(
                                                            newRegionLocation));
                                                    newRegion.setRegionLocation(newRegionLocation);
                                                    hearsayDAOBean.getRegionDAO().save(newRegion);

                                                    utrRegionList.add(newRegion);
                                                }

                                            }

                                            if (strandType.equals(StrandType.PLUS)) {

                                                if (transcriptRange.containsInteger(proteinLocation.getStart())) {

                                                    transcriptLocation.setStop(proteinLocation.getStop() - 1);
                                                    hearsayDAOBean.getLocationDAO().save(transcriptLocation);

                                                    regionLocation.setStop(regionStart + transcriptLocation.diff());
                                                    hearsayDAOBean.getLocationDAO().save(regionLocation);

                                                    Region newRegion = new Region(RegionType.EXON);
                                                    newRegion.setId(hearsayDAOBean.getRegionDAO().save(newRegion));
                                                    newRegion.setAlignment(alignment);

                                                    Location newTranscriptLocation = new Location(alignment
                                                            .getProteinLocation().getStart(), transcriptStop);
                                                    newTranscriptLocation.setId(hearsayDAOBean.getLocationDAO().save(
                                                            newTranscriptLocation));
                                                    newRegion.setTranscriptLocation(newTranscriptLocation);

                                                    Location newRegionLocation = new Location(
                                                            regionLocation.getStop() + 1, regionLocation.getStop() + 1
                                                                    + newTranscriptLocation.diff());
                                                    newRegionLocation.setId(hearsayDAOBean.getLocationDAO().save(
                                                            newRegionLocation));
                                                    newRegion.setRegionLocation(newRegionLocation);
                                                    hearsayDAOBean.getRegionDAO().save(newRegion);

                                                    utrRegionList.add(newRegion);

                                                }

                                                if (transcriptRange.containsInteger(proteinLocation.getStop())) {

                                                    transcriptLocation.setStop(proteinLocation.getStop());
                                                    hearsayDAOBean.getLocationDAO().save(transcriptLocation);

                                                    regionLocation.setStop(regionStart + transcriptLocation.diff());
                                                    hearsayDAOBean.getLocationDAO().save(regionLocation);

                                                    Region newRegion = new Region(RegionType.EXON);
                                                    newRegion.setId(hearsayDAOBean.getRegionDAO().save(newRegion));
                                                    newRegion.setAlignment(alignment);

                                                    Location newTranscriptLocation = new Location(
                                                            proteinLocation.getStop() + 1, transcriptStop);
                                                    newTranscriptLocation.setId(hearsayDAOBean.getLocationDAO().save(
                                                            newTranscriptLocation));
                                                    newRegion.setTranscriptLocation(newTranscriptLocation);

                                                    Location newRegionLocation = new Location(regionStop
                                                            - newTranscriptLocation.diff(), regionStop);
                                                    newRegionLocation.setId(hearsayDAOBean.getLocationDAO().save(
                                                            newRegionLocation));
                                                    newRegion.setRegionLocation(newRegionLocation);
                                                    hearsayDAOBean.getRegionDAO().save(newRegion);

                                                    utrRegionList.add(newRegion);

                                                }

                                            }

                                        }

                                        regionList.addAll(utrRegionList);

                                        Collections.sort(regionList, new Comparator<Region>() {
                                            @Override
                                            public int compare(Region o1, Region o2) {
                                                return Integer.compare(o1.getRegionLocation().getStart(), o2
                                                        .getRegionLocation().getStart());
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

                                            hearsayDAOBean.getRegionDAO().save(region);
                                        }

                                        // for (Region region : regionList) {
                                        // Location regionLocation = region.getRegionLocation();
                                        // Location transcriptLocation = region.getTranscriptLocation();
                                        // if (transcriptLocation == null) {
                                        // logger.info("{}, {}", region.toString(), regionLocation.toString());
                                        // continue;
                                        // }
                                        // logger.info("{}, {}, {}", region.toString(), regionLocation.toString(),
                                        // transcriptLocation.toString());
                                        // }

                                        List<Region> intronRegionList = new ArrayList<Region>();

                                        // adding intron regions
                                        Region previousRegion = null;
                                        for (Region currentRegion : regionList) {
                                            if (previousRegion == null) {
                                                previousRegion = currentRegion;
                                                continue;
                                            }

                                            if (previousRegion.getRegionLocation().getStop()
                                                    .equals(currentRegion.getRegionLocation().getStart() - 1)) {
                                                previousRegion = currentRegion;
                                                continue;
                                            }

                                            Region region = new Region(RegionType.INTRON);
                                            region.setId(hearsayDAOBean.getRegionDAO().save(region));
                                            region.setAlignment(alignment);

                                            Location regionLocation = new Location(previousRegion.getRegionLocation()
                                                    .getStop() + 1, currentRegion.getRegionLocation().getStart() - 1);
                                            regionLocation.setId(hearsayDAOBean.getLocationDAO().save(regionLocation));
                                            region.setRegionLocation(regionLocation);
                                            hearsayDAOBean.getRegionDAO().save(region);
                                            intronRegionList.add(region);

                                            previousRegion = currentRegion;
                                        }

                                        regionList.addAll(intronRegionList);

                                        Collections.sort(regionList, new Comparator<Region>() {
                                            @Override
                                            public int compare(Region o1, Region o2) {
                                                return Integer.compare(o1.getRegionLocation().getStart(), o2
                                                        .getRegionLocation().getStart());
                                            }
                                        });

                                        for (Region region : regionList) {
                                            Location regionLocation = region.getRegionLocation();
                                            Location transcriptLocation = region.getTranscriptLocation();
                                            if (transcriptLocation == null) {
                                                logger.debug("{}, {}", region.toString(), regionLocation.toString());
                                                continue;
                                            }
                                            logger.debug("{}, {}, {}", region.toString(), regionLocation.toString(),
                                                    transcriptLocation.toString());
                                        }

                                        alignment.setRegions(regionList);

                                        hearsayDAOBean.getAlignmentDAO().save(alignment);

                                    }

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

    public HearsayDAOBean getHearsayDAOBean() {
        return hearsayDAOBean;
    }

    public void setHearsayDAOBean(HearsayDAOBean hearsayDAOBean) {
        this.hearsayDAOBean = hearsayDAOBean;
    }

}
