package org.renci.hearsay.commands.ncbi;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.karaf.shell.api.action.Command;
import org.apache.karaf.shell.api.action.Action;
import org.renci.hearsay.dao.HearsayDAOBean;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

@Command(scope = "ncbi", name = "pull", description = "Pull Everything")
public class PullAction implements Action {

    private final Logger logger = LoggerFactory.getLogger(PullAction.class);

    private HearsayDAOBean hearsayDAOBean;

    public PullAction() {
        super();
    }

    @Override
    public Object execute() {
        logger.debug("ENTERING execute()");
        PullRunnable runnable = new PullRunnable();
        runnable.setHearsayDAOBean(hearsayDAOBean);
        ExecutorService es = Executors.newSingleThreadExecutor();
        es.submit(runnable);
        es.shutdown();
        return null;
    }

    public HearsayDAOBean getHearsayDAOBean() {
        return hearsayDAOBean;
    }

    public void setHearsayDAOBean(HearsayDAOBean hearsayDAOBean) {
        this.hearsayDAOBean = hearsayDAOBean;
    }

}
