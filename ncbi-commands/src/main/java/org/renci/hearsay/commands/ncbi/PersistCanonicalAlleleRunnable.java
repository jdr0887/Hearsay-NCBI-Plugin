package org.renci.hearsay.commands.ncbi;

import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.math.NumberUtils;
import org.renci.clinvar.MeasureSetType;
import org.renci.clinvar.MeasureSetType.Measure;
import org.renci.clinvar.MeasureSetType.Measure.AttributeSet;
import org.renci.clinvar.MeasureSetType.Measure.AttributeSet.Attribute;
import org.renci.clinvar.ReferenceAssertionType;
import org.renci.clinvar.SetElementSetType;
import org.renci.clinvar.XrefType;
import org.renci.hearsay.dao.HearsayDAOBean;
import org.renci.hearsay.dao.model.CanonicalAllele;
import org.renci.hearsay.dao.model.ComplexityType;
import org.renci.hearsay.dao.model.Identifier;
import org.renci.hearsay.dao.model.IntronOffset;
import org.renci.hearsay.dao.model.Location;
import org.renci.hearsay.dao.model.MolecularConsequence;
import org.renci.hearsay.dao.model.MolecularConsequenceType;
import org.renci.hearsay.dao.model.ReferenceCoordinate;
import org.renci.hearsay.dao.model.ReferenceSequence;
import org.renci.hearsay.dao.model.SimpleAllele;
import org.renci.hearsay.dao.model.SimpleAlleleType;
import org.renci.hearsay.dao.model.StrandType;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class PersistCanonicalAlleleRunnable implements Runnable {

    private static final Logger logger = LoggerFactory.getLogger(PersistCanonicalAlleleRunnable.class);

    private HearsayDAOBean hearsayDAOBean;

    private static final Pattern locationPattern = Pattern.compile("[A-Za-z]");

    private ReferenceAssertionType rat;

    public PersistCanonicalAlleleRunnable(HearsayDAOBean hearsayDAOBean, ReferenceAssertionType rat) {
        super();
        this.hearsayDAOBean = hearsayDAOBean;
        this.rat = rat;
    }

    @Override
    public void run() {

        try {
            ReferenceAssertionType.ClinVarAccession clinVarAccession = rat.getClinVarAccession();

            CanonicalAllele canonicalAllele = new CanonicalAllele();
            canonicalAllele.setActive("current".equals(rat.getRecordStatus()));
            canonicalAllele.setVersion(clinVarAccession.getVersion().toString());

            MeasureSetType mst = rat.getMeasureSet();
            List<SetElementSetType> measureSetTypeNameList = mst.getName();
            for (SetElementSetType mstNameList : measureSetTypeNameList) {
                if ("Preferred".equals(mstNameList.getElementValue().getType())) {
                    canonicalAllele.setName(mstNameList.getElementValue().getValue());
                    break;
                }
            }

            ComplexityType complexityType = ComplexityType.SIMPLE;
            if (mst.getMeasure().size() > 1) {
                complexityType = ComplexityType.COMPLEX;
            }
            canonicalAllele.setComplexityType(complexityType);

            canonicalAllele.setId(hearsayDAOBean.getCanonicalAlleleDAO().save(canonicalAllele));

            if ("Variant".equals(mst.getType())) {
                Identifier identifier = new Identifier("www.ncbi.nlm.nih.gov/clinvar/variation", mst.getID().toString());
                identifier.setId(hearsayDAOBean.getIdentifierDAO().save(identifier));
                logger.info(identifier.toString());
                canonicalAllele.getRelatedIdentifiers().add(identifier);
            }

            // for (MoleculeType mType : MoleculeType.values()) {
            // if (mType.getPrefixes().contains(canonicalAllele.getName().substring(0, 3))) {
            // canonicalAllele.setMoleculeType(mType);
            // break;
            // }
            // }

            for (Measure measure : mst.getMeasure()) {

                List<AttributeSet> attributeSetList = measure.getAttributeSet();
                
                String dbSNPId = null;
                for (XrefType xref : measure.getXRef()) {
                    if ("dbSNP".equals(xref.getDB())) {
                        dbSNPId = xref.getID();
                        break;
                    }
                }

                Set<SimpleAllele> simpleAlleleSet = new HashSet<SimpleAllele>();
                for (AttributeSet attributeSet : attributeSetList) {
                    Attribute attribute = attributeSet.getAttribute();
                    for (SimpleAlleleType saType : SimpleAlleleType.values()) {
                        if (!saType.getType().equals(attribute.getType())) {
                            continue;
                        }
                        SimpleAllele simpleAllele = new SimpleAllele();
                        simpleAllele.setName(attribute.getValue());
                        if (saType.equals(SimpleAlleleType.TRANSCRIPT) || saType.equals(SimpleAlleleType.GENOMIC)) {
                            simpleAllele.setAllele(attribute.getValue().substring(attribute.getValue().length() - 1,
                                    attribute.getValue().length()));
                        }
                        simpleAllele.setType(saType);
                        simpleAllele.setCanonicalAllele(canonicalAllele);
                        simpleAllele.setId(hearsayDAOBean.getSimpleAlleleDAO().save(simpleAllele));
                        simpleAlleleSet.add(simpleAllele);
                        break;
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
                                if (sa.getName().startsWith(refSeqId.substring(0, refSeqId.indexOf(":")))) {
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
                    String referenceSequenceAccession = sa.getName().substring(0, sa.getName().indexOf(":"));
                    ReferenceSequence referenceSequence = null;
                    List<ReferenceSequence> potentialRefSeqList = hearsayDAOBean.getReferenceSequenceDAO()
                            .findByIdentifierValue(referenceSequenceAccession);
                    if (CollectionUtils.isNotEmpty(potentialRefSeqList)) {
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
                                    intron.setId(hearsayDAOBean.getIntronOffsetDAO().save(intron));
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
                                    intron.setId(hearsayDAOBean.getIntronOffsetDAO().save(intron));
                                    referenceCoordinate.setIntronOffset(intron);
                                }

                            }
                        }
                    }
                    referenceCoordinate.setId(hearsayDAOBean.getReferenceCoordinateDAO().save(referenceCoordinate));

                    if (StringUtils.isNotEmpty(dbSNPId)) {
                        Identifier identifier = new Identifier("www.ncbi.nlm.nih.gov/snp", dbSNPId);
                        List<Identifier> possibleIdentifiers = hearsayDAOBean.getIdentifierDAO().findByExample(
                                identifier);
                        if (CollectionUtils.isNotEmpty(possibleIdentifiers)) {
                            identifier = possibleIdentifiers.get(0);
                        } else {
                            identifier.setId(hearsayDAOBean.getIdentifierDAO().save(identifier));
                        }
                        logger.info(identifier.toString());
                        referenceCoordinate.getIdentifiers().add(identifier);
                    }
                    hearsayDAOBean.getReferenceCoordinateDAO().save(referenceCoordinate);

                    sa.setReferenceCoordinate(referenceCoordinate);
                    hearsayDAOBean.getSimpleAlleleDAO().save(sa);
                }

                canonicalAllele.getRelatedSimpleAlleles().addAll(simpleAlleleSet);

            }
            hearsayDAOBean.getCanonicalAlleleDAO().save(canonicalAllele);

        } catch (Exception e) {
            logger.error("Error", e);
            e.printStackTrace();
        }

    }

}
