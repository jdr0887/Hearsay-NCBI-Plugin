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

@Command(scope = "ncbi", name = "pull-genome-references", description = "Pull Genome References")
@Service
public class PullGenomeReferencesAction implements Action {

    private final Logger logger = LoggerFactory.getLogger(PullGenomeReferencesAction.class);

    @Reference
    private HearsayDAOBeanService hearsayDAOBeanService;

    public PullGenomeReferencesAction() {
        super();
    }

    @Override
    public Object execute() {
        logger.debug("ENTERING execute()");
        ExecutorService es = Executors.newSingleThreadExecutor();
        es.submit(new PullGenomeReferencesRunnable(hearsayDAOBeanService));
        es.shutdown();
        return null;
    }

}
