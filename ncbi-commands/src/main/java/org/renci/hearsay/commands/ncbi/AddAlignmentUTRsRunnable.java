package org.renci.hearsay.commands.ncbi;

import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.collections.CollectionUtils;
import org.renci.hearsay.dao.HearsayDAOBean;
import org.renci.hearsay.dao.model.Alignment;
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

            ExecutorService es = Executors.newFixedThreadPool(4);

            List<Alignment> alignmentList = hearsayDAOBean.getAlignmentDAO().findAll();

            if (CollectionUtils.isNotEmpty(alignmentList)) {

                logger.debug("alignmentList.size(): {}", alignmentList.size());

                for (Alignment alignment : alignmentList) {

                    if (CollectionUtils.isNotEmpty(alignment.getRegions())) {
                        es.submit(new AddUTRRunnable(hearsayDAOBean, alignment));
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

    public HearsayDAOBean getHearsayDAOBean() {
        return hearsayDAOBean;
    }

    public void setHearsayDAOBean(HearsayDAOBean hearsayDAOBean) {
        this.hearsayDAOBean = hearsayDAOBean;
    }

}
