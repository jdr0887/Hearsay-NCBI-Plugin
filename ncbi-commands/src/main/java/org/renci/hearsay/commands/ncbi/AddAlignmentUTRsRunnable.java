package org.renci.hearsay.commands.ncbi;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

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
            if (geneList != null && !geneList.isEmpty()) {
                for (Gene gene : geneList) {

                    logger.info(gene.toString());

                    List<ReferenceSequence> referenceSequenceList = hearsayDAOBean.getReferenceSequenceDAO()
                            .findByGeneId(gene.getId());

                    if (referenceSequenceList != null && !referenceSequenceList.isEmpty()) {
                        
                        for (ReferenceSequence referenceSequence : referenceSequenceList) {
                            
                            StrandType strandType = referenceSequence.getStrandType();

                            logger.info(referenceSequence.toString());

                            List<Alignment> alignmentList = hearsayDAOBean.getAlignmentDAO().findByReferenceSequenceId(
                                    referenceSequence.getId());

                            if (alignmentList != null && !alignmentList.isEmpty()) {

                                logger.info("alignmentList.size(): {}", alignmentList.size());

                                for (Alignment alignment : alignmentList) {

                                    // adding utr regions
                                    List<Region> utrRegionList = new ArrayList<Region>();
                                    Iterator<Region> regionIter = alignment.getRegions().iterator();
                                    Region firstRegion = alignment.getRegions().get(0);
                                    Region lastRegion = alignment.getRegions().get(alignment.getRegions().size() - 1);
                                    while (regionIter.hasNext()) {
                                        Region region = regionIter.next();

                                        if (alignment.getProteinLocation().getStart()
                                                .equals(firstRegion.getTranscriptLocation().getStart())
                                                && alignment.getProteinLocation().getStop()
                                                        .equals(lastRegion.getTranscriptLocation().getStop())) {
                                            continue;
                                        }

                                        Location regionLocation = region.getRegionLocation();
                                        int regionStart = regionLocation.getStart();
                                        int regionStop = regionLocation.getStop();

                                        Location transcriptLocation = region.getTranscriptLocation();
                                        IntRange transcriptRange = transcriptLocation.toRange();

                                        int transcriptStart = transcriptLocation.getStart();
                                        int transcriptStop = transcriptLocation.getStop();

                                        if (strandType.equals(StrandType.MINUS)
                                                && transcriptRange.containsInteger(alignment.getProteinLocation()
                                                        .getStart())) {

                                            transcriptLocation.setStop(alignment.getProteinLocation().getStart());
                                            hearsayDAOBean.getLocationDAO().save(transcriptLocation);

                                            int diff = transcriptStart - transcriptLocation.getStop();
                                            regionLocation.setStop(regionLocation.getStart() + diff);
                                            hearsayDAOBean.getLocationDAO().save(regionLocation);

                                            Region newRegion = new Region(RegionType.UTR5);
                                            newRegion.setAlignment(alignment);

                                            Location newTranscriptLocation = new Location(alignment
                                                    .getProteinLocation().getStart() - 1, transcriptStop);
                                            newTranscriptLocation.setId(hearsayDAOBean.getLocationDAO().save(
                                                    newTranscriptLocation));
                                            newRegion.setTranscriptLocation(newTranscriptLocation);

                                            Location newRegionLocation = new Location(regionLocation.getStop() + 1,
                                                    regionLocation.getStop() + 1
                                                            + newRegion.getTranscriptLocation().diff());
                                            newRegionLocation.setId(hearsayDAOBean.getLocationDAO().save(
                                                    newRegionLocation));
                                            newRegion.setRegionLocation(newRegionLocation);

                                            newRegion.setId(hearsayDAOBean.getRegionDAO().save(newRegion));
                                            utrRegionList.add(newRegion);
                                        }

                                        if (strandType.equals(StrandType.MINUS)
                                                && transcriptRange.containsInteger(alignment.getProteinLocation()
                                                        .getStop())) {

                                            transcriptLocation.setStart(alignment.getProteinLocation().getStop());
                                            hearsayDAOBean.getLocationDAO().save(transcriptLocation);

                                            int diff = alignment.getProteinLocation().getStop() - transcriptStop;
                                            regionLocation.setStart(regionLocation.getStop() - diff);
                                            hearsayDAOBean.getLocationDAO().save(regionLocation);

                                            Region newRegion = new Region(RegionType.UTR3);
                                            newRegion.setAlignment(alignment);

                                            Location newTranscriptLocation = new Location(transcriptStart, alignment
                                                    .getProteinLocation().getStop() + 1);
                                            newTranscriptLocation.setId(hearsayDAOBean.getLocationDAO().save(
                                                    newTranscriptLocation));
                                            newRegion.setTranscriptLocation(newTranscriptLocation);

                                            Location newRegionLocation = new Location(regionLocation.getStart() - 2
                                                    - newTranscriptLocation.diff(), regionLocation.getStart() - 1);
                                            newRegionLocation.setId(hearsayDAOBean.getLocationDAO().save(
                                                    newRegionLocation));
                                            newRegion.setRegionLocation(newRegionLocation);

                                            newRegion.setId(hearsayDAOBean.getRegionDAO().save(newRegion));
                                            utrRegionList.add(newRegion);
                                        }

                                        if (strandType.equals(StrandType.PLUS)
                                                && transcriptRange.containsInteger(alignment.getProteinLocation()
                                                        .getStart())) {

                                            transcriptLocation.setStart(alignment.getProteinLocation().getStart());
                                            hearsayDAOBean.getLocationDAO().save(transcriptLocation);

                                            regionLocation.setStart(regionStop - transcriptLocation.diff());
                                            hearsayDAOBean.getLocationDAO().save(regionLocation);

                                            Region newRegion = new Region(RegionType.UTR5);
                                            newRegion.setAlignment(alignment);

                                            Location newTranscriptLocation = new Location(transcriptStart, alignment
                                                    .getProteinLocation().getStart() - 1);
                                            newTranscriptLocation.setId(hearsayDAOBean.getLocationDAO().save(
                                                    newTranscriptLocation));
                                            newRegion.setTranscriptLocation(newTranscriptLocation);

                                            Location newRegionLocation = new Location(regionLocation.getStart()
                                                    - newTranscriptLocation.diff(), regionLocation.getStart() - 1);
                                            newRegionLocation.setId(hearsayDAOBean.getLocationDAO().save(
                                                    newRegionLocation));
                                            newRegion.setRegionLocation(newRegionLocation);

                                            newRegion.setId(hearsayDAOBean.getRegionDAO().save(newRegion));
                                            utrRegionList.add(newRegion);

                                        }

                                        if (strandType.equals(StrandType.PLUS)
                                                && transcriptRange.containsInteger(alignment.getProteinLocation()
                                                        .getStop())) {

                                            transcriptLocation.setStop(alignment.getProteinLocation().getStop());
                                            hearsayDAOBean.getLocationDAO().save(transcriptLocation);

                                            regionLocation.setStop(regionStart + transcriptLocation.diff());
                                            hearsayDAOBean.getLocationDAO().save(regionLocation);

                                            Region newRegion = new Region(RegionType.UTR3);
                                            newRegion.setAlignment(alignment);

                                            Location newTranscriptLocation = new Location(alignment
                                                    .getProteinLocation().getStop() + 1, transcriptStop);
                                            newTranscriptLocation.setId(hearsayDAOBean.getLocationDAO().save(
                                                    newTranscriptLocation));
                                            newRegion.setTranscriptLocation(newTranscriptLocation);

                                            Location newRegionLocation = new Location(regionStop
                                                    - newTranscriptLocation.diff(), regionStop);
                                            newRegionLocation.setId(hearsayDAOBean.getLocationDAO().save(
                                                    newRegionLocation));
                                            newRegion.setRegionLocation(newRegionLocation);

                                            newRegion.setId(hearsayDAOBean.getRegionDAO().save(newRegion));
                                            utrRegionList.add(newRegion);

                                        }

                                    }

                                    alignment.getRegions().addAll(utrRegionList);
                                    hearsayDAOBean.getAlignmentDAO().save(alignment);

                                    for (Region region : alignment.getRegions()) {

                                        if (strandType.equals(StrandType.PLUS)
                                                && region.getTranscriptLocation().getStop() < alignment
                                                        .getProteinLocation().getStart()) {
                                            region.setRegionType(RegionType.UTR5);
                                        }

                                        if (strandType.equals(StrandType.PLUS)
                                                && region.getTranscriptLocation().getStop() > alignment
                                                        .getProteinLocation().getStop()) {
                                            region.setRegionType(RegionType.UTR3);
                                        }

                                        if (strandType.equals(StrandType.MINUS)
                                                && region.getTranscriptLocation().getStop() < alignment
                                                        .getProteinLocation().getStart()) {
                                            region.setRegionType(RegionType.UTR5);
                                        }

                                        if (strandType.equals(StrandType.MINUS)
                                                && region.getTranscriptLocation().getStop() > alignment
                                                        .getProteinLocation().getStop()) {
                                            region.setRegionType(RegionType.UTR3);
                                        }

                                        hearsayDAOBean.getRegionDAO().save(region);
                                    }

                                    logger.info("Adding INTRONs");

                                    // adding intron regions
                                    List<Region> intronRegionList = new ArrayList<Region>();
                                    regionIter = alignment.getRegions().iterator();
                                    Region previousRegion = null;
                                    while (regionIter.hasNext()) {
                                        Region currentRegion = regionIter.next();
                                        if (previousRegion == null) {
                                            previousRegion = currentRegion;
                                            continue;
                                        }

                                        if (!currentRegion.getRegionLocation().getStart()
                                                .equals(previousRegion.getRegionLocation().getStop() - 1)) {
                                            Region region = new Region(RegionType.INTRON);
                                            region.setAlignment(alignment);
                                            Location regionLocation = new Location(previousRegion.getRegionLocation()
                                                    .getStop() + 1, currentRegion.getRegionLocation().getStart() - 1);
                                            regionLocation.setId(hearsayDAOBean.getLocationDAO().save(regionLocation));
                                            region.setRegionLocation(regionLocation);
                                            region.setId(hearsayDAOBean.getRegionDAO().save(region));
                                            intronRegionList.add(region);
                                        }

                                        previousRegion = currentRegion;
                                    }

                                    alignment.getRegions().addAll(intronRegionList);
                                    hearsayDAOBean.getAlignmentDAO().save(alignment);

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
