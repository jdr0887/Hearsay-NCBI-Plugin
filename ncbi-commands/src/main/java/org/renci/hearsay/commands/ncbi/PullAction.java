package org.renci.hearsay.commands.ncbi;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.karaf.shell.commands.Command;
import org.apache.karaf.shell.console.AbstractAction;
import org.renci.hearsay.dao.HearsayDAOBean;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

@Command(scope = "ncbi", name = "pull", description = "Pull Everything")
public class PullAction extends AbstractAction {

    private final Logger logger = LoggerFactory.getLogger(PullAction.class);

    private HearsayDAOBean hearsayDAOBean;

    public PullAction() {
        super();
    }

    @Override
    public Object doExecute() {
        logger.debug("ENTERING doExecute()");
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