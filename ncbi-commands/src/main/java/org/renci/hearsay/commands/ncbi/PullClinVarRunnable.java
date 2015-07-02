package org.renci.hearsay.commands.ncbi;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.math.NumberUtils;
import org.renci.clinvar.MeasureSetType;
import org.renci.clinvar.MeasureSetType.Measure;
import org.renci.clinvar.MeasureSetType.Measure.AttributeSet;
import org.renci.clinvar.MeasureSetType.Measure.AttributeSet.Attribute;
import org.renci.clinvar.PublicSetType;
import org.renci.clinvar.ReferenceAssertionType;
import org.renci.clinvar.ReleaseType;
import org.renci.clinvar.SetElementSetType;
import org.renci.clinvar.XrefType;
import org.renci.hearsay.commands.ncbi.util.FTPUtil;
import org.renci.hearsay.dao.HearsayDAOBean;
import org.renci.hearsay.dao.HearsayDAOException;
import org.renci.hearsay.dao.model.CanonicalAllele;
import org.renci.hearsay.dao.model.ComplexityType;
import org.renci.hearsay.dao.model.Identifier;
import org.renci.hearsay.dao.model.IntronOffset;
import org.renci.hearsay.dao.model.Location;
import org.renci.hearsay.dao.model.MolecularConsequence;
import org.renci.hearsay.dao.model.MolecularConsequenceType;
import org.renci.hearsay.dao.model.MoleculeType;
import org.renci.hearsay.dao.model.ReferenceCoordinate;
import org.renci.hearsay.dao.model.ReferenceSequence;
import org.renci.hearsay.dao.model.SimpleAllele;
import org.renci.hearsay.dao.model.SimpleAlleleType;
import org.renci.hearsay.dao.model.StrandType;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class PullClinVarRunnable implements Runnable {

    private final Logger logger = LoggerFactory.getLogger(PullClinVarRunnable.class);

    private HearsayDAOBean hearsayDAOBean;

    private static final Pattern locationPattern = Pattern.compile("[A-Za-z]");

    public PullClinVarRunnable() {
        super();
    }

    @Override
    public void run() {

        File clinvarDownload = FTPUtil.ncbiDownload("/pub/clinvar/xml", "ClinVarFullRelease_00-latest.xml.gz");

        try {
            JAXBContext jc = JAXBContext.newInstance(ReleaseType.class);
            Unmarshaller u = jc.createUnmarshaller();
            ReleaseType releaseType = (ReleaseType) u.unmarshal(new GZIPInputStream(
                    new FileInputStream(clinvarDownload)));
            List<PublicSetType> publicSetType = releaseType.getClinVarSet();
            for (PublicSetType pst : publicSetType) {
                ReferenceAssertionType rat = pst.getReferenceClinVarAssertion();
                ReferenceAssertionType.ClinVarAccession clinVarAccession = rat.getClinVarAccession();

                CanonicalAllele canonicalAllele = new CanonicalAllele();
                canonicalAllele.setActive("current".equals(rat.getRecordStatus()));
                canonicalAllele.setVersion(clinVarAccession.getVersion().toString());

                MeasureSetType mst = rat.getMeasureSet();
                List<SetElementSetType> measureSetTypeNameList = mst.getName();
                for (SetElementSetType mstNameList : measureSetTypeNameList) {
                    if ("Preferred".equals(mstNameList.getElementValue().getType())) {
                        canonicalAllele.setName(mstNameList.getElementValue().getValue());
                    }
                }

                ComplexityType complexityType = ComplexityType.SIMPLE;
                if (mst.getMeasure().size() > 1) {
                    complexityType = ComplexityType.COMPLEX;
                }
                canonicalAllele.setComplexityType(complexityType);

                canonicalAllele.setId(hearsayDAOBean.getCanonicalAlleleDAO().save(canonicalAllele));

                if ("Variant".equals(mst.getType())) {
                    Identifier identifier = new Identifier();
                    identifier.setSystem("http://www.ncbi.nlm.nih.gov/clinvar");
                    identifier.setValue(mst.getID().toString());
                    identifier.setId(hearsayDAOBean.getIdentifierDAO().save(identifier));
                    logger.info(identifier.toString());
                    canonicalAllele.getRelatedIdentifiers().add(identifier);
                }

                for (MoleculeType mType : MoleculeType.values()) {
                    if (mType.getPrefixes().contains(canonicalAllele.getName().substring(0, 3))) {
                        canonicalAllele.setMoleculeType(mType);
                        break;
                    }
                }

                for (Measure measure : mst.getMeasure()) {

                    List<AttributeSet> attributeSetList = measure.getAttributeSet();
                    String dbSNPId = null;
                    for (XrefType xref : measure.getXRef()) {
                        if ("dbSNP".equals(xref.getDB())) {
                            dbSNPId = xref.getID();
                        }
                    }

                    Set<SimpleAllele> simpleAlleleSet = new HashSet<SimpleAllele>();
                    for (AttributeSet attributeSet : attributeSetList) {
                        Attribute attribute = attributeSet.getAttribute();
                        for (SimpleAlleleType saType : SimpleAlleleType.values()) {
                            if (saType.getType().equals(attribute.getType())) {
                                SimpleAllele simpleAllele = new SimpleAllele();
                                simpleAllele.setName(attribute.getValue());
                                if (saType.equals(SimpleAlleleType.TRANSCRIPT)
                                        || saType.equals(SimpleAlleleType.GENOMIC)) {
                                    simpleAllele.setAllele(attribute.getValue().substring(
                                            attribute.getValue().length() - 1, attribute.getValue().length()));
                                }
                                simpleAllele.setType(saType);
                                simpleAlleleSet.add(simpleAllele);
                                break;
                            }
                        }
                    }

                    for (AttributeSet attributeSet : attributeSetList) {
                        Attribute attribute = attributeSet.getAttribute();
                        if ("MolecularConsequence".equals(attribute.getType())) {
                            String sequenceOntologyId = null;
                            String refSeqId = null;
                            for (XrefType xref : attributeSet.getXRef()) {
                                if ("Sequence Ontology".equals(xref.getDB())) {
                                    sequenceOntologyId = xref.getID();
                                }
                                if ("RefSeq".equals(xref.getDB())) {
                                    refSeqId = xref.getID();
                                }
                            }

                            if (StringUtils.isNotEmpty(refSeqId) && StringUtils.isNotBlank(sequenceOntologyId)) {
                                for (SimpleAllele sa : simpleAlleleSet) {
                                    if (sa.getName().startsWith(refSeqId.substring(0, refSeqId.indexOf(".")))) {
                                        MolecularConsequence mc = new MolecularConsequence(
                                                Integer.valueOf(sequenceOntologyId.replace("SO:", "")),
                                                MolecularConsequenceType.PRIMARY);
                                        mc.setId(hearsayDAOBean.getMolecularConsequenceDAO().save(mc));
                                        sa.getMolecularConsequences().add(mc);
                                    }
                                }
                            }
                        }
                    }

                    for (SimpleAllele sa : simpleAlleleSet) {
                        String referenceSequenceAccession = sa.getName().substring(0, sa.getName().indexOf("."));
                        ReferenceSequence referenceSequence = null;
                        List<ReferenceSequence> potentialRefSeqList = hearsayDAOBean.getReferenceSequenceDAO()
                                .findByIdentifierValue(referenceSequenceAccession);
                        if (potentialRefSeqList != null && !potentialRefSeqList.isEmpty()) {
                            referenceSequence = potentialRefSeqList.get(0);
                        }

                        ReferenceCoordinate referenceCoordinate = new ReferenceCoordinate();
                        referenceCoordinate.setReferenceSequence(referenceSequence);

                        String hgvsDescription = sa.getName().substring(sa.getName().indexOf(":") + 1,
                                sa.getName().length());
                        String type = hgvsDescription.substring(0, 1);
                        if (type.equals("c") || type.equals("g")) {
                            String s = hgvsDescription.substring(2);
                            Matcher locationMatcher = locationPattern.matcher(s);
                            locationMatcher.find();
                            String location = s.substring(0, locationMatcher.start());
                            s = s.substring(locationMatcher.start());
                            if (s.contains(">")) {
                                referenceCoordinate.setRefAllele(s.substring(0, s.indexOf(">")));
                                if (NumberUtils.isNumber(location)) {
                                    // a change in the coding
                                    Location rcLocation = new Location(Integer.valueOf(location) - 1,
                                            Integer.valueOf(location));
                                    rcLocation.setId(hearsayDAOBean.getLocationDAO().save(rcLocation));
                                    referenceCoordinate.setLocation(rcLocation);
                                } else {
                                    if (location.contains("-") && !location.startsWith("-")) {
                                        // a change in the 3' end of an intron
                                        Integer start = Integer.valueOf(location.substring(0, location.indexOf("-")));
                                        Integer end = Integer.valueOf(location.substring(location.indexOf("-") + 1,
                                                location.length()));
                                        IntronOffset intron = new IntronOffset(start, end, StrandType.MINUS);
                                        referenceCoordinate.setIntronOffset(intron);
                                    }

                                    if (location.startsWith("-")) {
                                        // a change 5' of the ATG (in the 5'UTR)
                                    }

                                    if (location.startsWith("*")) {
                                        // a change 3' of the stop codon (in the 3'UTR)
                                    }

                                    if (location.contains("+")) {
                                        // a change in the 5' end of an intron
                                        Integer end = Integer.valueOf(location.substring(location.indexOf("+") + 1,
                                                location.length()));
                                        Integer start = end - 1;
                                        IntronOffset intron = new IntronOffset(start, end, StrandType.PLUS);
                                        referenceCoordinate.setIntronOffset(intron);
                                    }

                                }
                            }
                        }

                        if (StringUtils.isNotEmpty(dbSNPId)) {
                            Identifier identifier = new Identifier("www.ncbi.nlm.nih.gov/snp", dbSNPId);
                            List<Identifier> possibleIdentifiers = hearsayDAOBean.getIdentifierDAO().findByExample(
                                    identifier);
                            if (possibleIdentifiers != null && !possibleIdentifiers.isEmpty()) {
                                identifier = possibleIdentifiers.get(0);
                            } else {
                                identifier.setId(hearsayDAOBean.getIdentifierDAO().save(identifier));
                            }
                            logger.info(identifier.toString());
                            referenceCoordinate.getIdentifiers().add(identifier);
                        }

                        sa.setReferenceCoordinate(referenceCoordinate);
                        hearsayDAOBean.getSimpleAlleleDAO().save(sa);
                    }

                    canonicalAllele.getRelatedSimpleAlleles().addAll(simpleAlleleSet);

                }

            }

        } catch (HearsayDAOException | JAXBException | IOException e) {
            logger.error("JAXBException | IOException", e);
        }

    }

    public HearsayDAOBean getHearsayDAOBean() {
        return hearsayDAOBean;
    }

    public void setHearsayDAOBean(HearsayDAOBean hearsayDAOBean) {
        this.hearsayDAOBean = hearsayDAOBean;
    }

}
