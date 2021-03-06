package org.renci.hearsay.commands.ncbi;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.karaf.shell.api.action.Action;
import org.apache.karaf.shell.api.action.Command;
import org.apache.karaf.shell.api.action.lifecycle.Reference;
import org.apache.karaf.shell.api.action.lifecycle.Service;
import org.renci.hearsay.dao.HearsayDAOBeanService;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

@Command(scope = "ncbi", name = "pull-genes", description = "Pull Genes")
@Service
public class PullGenesAction implements Action {

    private final Logger logger = LoggerFactory.getLogger(PullGenesAction.class);

    @Reference
    private HearsayDAOBeanService hearsayDAOBeanService;

    public PullGenesAction() {
        super();
    }

    @Override
    public Object execute() {
        logger.debug("ENTERING execute()");
        ExecutorService es = Executors.newSingleThreadExecutor();
        es.submit(new PullGenesRunnable(hearsayDAOBeanService));
        es.shutdown();
        return null;
    }

}
