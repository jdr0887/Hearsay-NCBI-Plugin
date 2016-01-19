package org.renci.hearsay.commands.ncbi;

import static org.renci.hearsay.commands.ncbi.Constants.IDENTIFIER_KEY_NUCCORE;

import java.util.Arrays;
import java.util.Formatter;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Set;

import org.apache.commons.collections4.CollectionUtils;
import org.apache.karaf.shell.api.action.Action;
import org.apache.karaf.shell.api.action.Command;
import org.apache.karaf.shell.api.action.lifecycle.Reference;
import org.apache.karaf.shell.api.action.lifecycle.Service;
import org.renci.hearsay.dao.HearsayDAOException;
import org.renci.hearsay.dao.IdentifierDAO;
import org.renci.hearsay.dao.ReferenceSequenceDAO;
import org.renci.hearsay.dao.model.Identifier;
import org.renci.hearsay.dao.model.ReferenceSequence;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

@Command(scope = "ncbi", name = "list-reference-sequences", description = "List Reference Sequences")
@Service
public class ListReferenceSequencesAction implements Action {

    private final Logger logger = LoggerFactory.getLogger(ListReferenceSequencesAction.class);

    @Reference
    private IdentifierDAO identifierDAO;

    @Reference
    private ReferenceSequenceDAO referenceSequenceDAO;

    public ListReferenceSequencesAction() {
        super();
    }

    @Override
    public Object execute() {
        logger.debug("ENTERING execute()");

        try {
            String format = "%1$-10s %2$-20s %3$-20s %4$-20s %5$-12s %6$s%n";
            StringBuilder sb = new StringBuilder();
            Formatter formatter = new Formatter(sb, Locale.US);
            formatter.format(format, "ID", "Accession", "Type");

            List<Identifier> identiferList = identifierDAO.findByExample(new Identifier(IDENTIFIER_KEY_NUCCORE, null));
            Set<ReferenceSequence> referenceSequences = new HashSet<ReferenceSequence>();
            if (CollectionUtils.isNotEmpty(identiferList)) {
                for (Identifier identifier : identiferList) {
                    List<ReferenceSequence> refSeqs = referenceSequenceDAO.findByIdentifiers(Arrays.asList(identifier));
                    if (CollectionUtils.isNotEmpty(refSeqs)) {
                        referenceSequences.addAll(refSeqs);
                    }
                }
            }

            if (CollectionUtils.isNotEmpty(referenceSequences)) {

                for (ReferenceSequence refSeq : referenceSequences) {

                    if (CollectionUtils.isNotEmpty(refSeq.getIdentifiers())) {
                        String rnaAccession = "";
                        String proteinAccession = "";
                        String genomicAccession = "";
                        for (Identifier identifier : refSeq.getIdentifiers()) {
                            if (identifier.getSystem().endsWith("nuccore")) {
                                rnaAccession = identifier.getValue();
                            }
                            if (identifier.getSystem().endsWith("protein")) {
                                proteinAccession = identifier.getValue();
                            }
                            if (identifier.getSystem().endsWith("genome")) {
                                genomicAccession = identifier.getValue();
                            }
                        }
                        formatter.format(format, refSeq.getId().toString(), rnaAccession, proteinAccession, genomicAccession,
                                refSeq.getStrandType().toString(), refSeq.getType().toString());
                    }

                    formatter.flush();
                }

            }
            System.out.println(formatter.toString());
            formatter.close();
        } catch (HearsayDAOException e) {
            e.printStackTrace();
        }

        return null;
    }

}
