package org.renci.hearsay.commands.ncbi;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.karaf.shell.api.action.Command;
import org.apache.karaf.shell.api.action.Action;
import org.renci.hearsay.dao.HearsayDAOBean;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

@Command(scope = "ncbi", name = "pull-genome-references", description = "Pull Genome References")
public class PullGenomeReferencesAction implements Action {

    private final Logger logger = LoggerFactory.getLogger(PullGenomeReferencesAction.class);

    private HearsayDAOBean hearsayDAOBean;

    public PullGenomeReferencesAction() {
        super();
    }

    @Override
    public Object execute() {
        logger.debug("ENTERING execute()");
        PullGenomeReferencesRunnable runnable = new PullGenomeReferencesRunnable();
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
