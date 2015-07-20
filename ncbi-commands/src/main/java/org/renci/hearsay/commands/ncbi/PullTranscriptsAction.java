package org.renci.hearsay.commands.ncbi;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.karaf.shell.commands.Command;
import org.apache.karaf.shell.console.AbstractAction;
import org.renci.hearsay.dao.HearsayDAOBean;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

@Command(scope = "ncbi", name = "pull-transcripts", description = "Pull Transcripts")
public class PullTranscriptsAction extends AbstractAction {

    private final Logger logger = LoggerFactory.getLogger(PullTranscriptsAction.class);

    private HearsayDAOBean hearsayDAOBean;

    public PullTranscriptsAction() {
        super();
    }

    @Override
    public Object doExecute() {
        logger.debug("ENTERING doExecute()");
        PullTranscriptsRunnable runnable = new PullTranscriptsRunnable();
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
