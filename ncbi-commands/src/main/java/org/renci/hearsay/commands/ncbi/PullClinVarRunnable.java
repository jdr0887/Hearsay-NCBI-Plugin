package org.renci.hearsay.commands.ncbi;

import static org.renci.hearsay.commands.ncbi.Constants.IDENTIFIER_KEY_NUCCORE;
import static org.renci.hearsay.commands.ncbi.Constants.IDENTIFIER_KEY_SNP;
import static org.renci.hearsay.commands.ncbi.Constants.IDENTIFIER_KEY_VARIATION;

import java.io.File;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Unmarshaller;

import org.apache.commons.collections4.CollectionUtils;
import org.renci.clinvar.MeasureSetType;
import org.renci.clinvar.MeasureSetType.Measure;
import org.renci.clinvar.MeasureSetType.Measure.AttributeSet;
import org.renci.clinvar.MeasureSetType.Measure.AttributeSet.Attribute;
import org.renci.clinvar.PublicSetType;
import org.renci.clinvar.ReferenceAssertionType;
import org.renci.clinvar.ReleaseType;
import org.renci.clinvar.SetElementSetType;
import org.renci.clinvar.SetElementSetType.ElementValue;
import org.renci.clinvar.XrefType;
import org.renci.hearsay.commands.ncbi.util.FTPUtil;
import org.renci.hearsay.dao.HearsayDAOBeanService;
import org.renci.hearsay.dao.model.CanonicalAllele;
import org.renci.hearsay.dao.model.CanonicalAlleleType;
import org.renci.hearsay.dao.model.ComplexityType;
import org.renci.hearsay.dao.model.ContextualAllele;
import org.renci.hearsay.dao.model.ContextualAlleleName;
import org.renci.hearsay.dao.model.ContextualAlleleNameType;
import org.renci.hearsay.dao.model.ContextualAlleleType;
import org.renci.hearsay.dao.model.DirectionType;
import org.renci.hearsay.dao.model.ExternalOffsetPosition;
import org.renci.hearsay.dao.model.Identifier;
import org.renci.hearsay.dao.model.ReferenceCoordinate;
import org.renci.hearsay.dao.model.ReferenceSequence;
import org.renci.hearsay.dao.model.ReferenceSequenceType;
import org.renci.hgvs.HGVSParser;
import org.renci.hgvs.model.AlleleInfo;
import org.renci.hgvs.model.dna.DNAChangeType;
import org.renci.hgvs.model.dna.DNAVariantMutation;
import org.renci.hgvs.model.dna.SubstitutionAlleleInfo;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class PullClinVarRunnable implements Runnable {

    private static final Logger logger = LoggerFactory.getLogger(PullClinVarRunnable.class);

    private static final List<String> allowedTypes = Arrays.asList("single nucleotide variant", "Duplication", "Deletion", "Indel",
            "inversion");

    private HearsayDAOBeanService hearsayDAOBeanService;

    public PullClinVarRunnable(HearsayDAOBeanService hearsayDAOBeanService) {
        super();
        this.hearsayDAOBeanService = hearsayDAOBeanService;
    }

    @Override
    public void run() {
        try {
            File clinvarDownload = FTPUtil.ncbiDownload("/pub/clinvar/xml", "ClinVarFullRelease_00-latest.xml.gz");
            JAXBContext jc = JAXBContext.newInstance(ReleaseType.class);
            Unmarshaller u = jc.createUnmarshaller();
            ReleaseType releaseType = (ReleaseType) u.unmarshal(new GZIPInputStream(new FileInputStream(clinvarDownload)));
            List<PublicSetType> publicSetTypeList = releaseType.getClinVarSet();

            if (CollectionUtils.isEmpty(publicSetTypeList)) {
                logger.warn("No PublicSetTypes found");
                return;
            }

            persistIdentifiers(publicSetTypeList);
            persistCanonicalAlleles(publicSetTypeList);
            // persistMeasureNameContextualAlleles(publicSetTypeList);
            persistMeasureAttributeContextualAlleles(publicSetTypeList);

        } catch (Exception e) {
            logger.error("Error", e);
        }

    }

    private void persistIdentifiers(List<PublicSetType> publicSetTypeList) {

        // has to be single threaded to avoid race condition
        for (PublicSetType pst : publicSetTypeList) {

            try {
                ReferenceAssertionType rat = pst.getReferenceClinVarAssertion();
                MeasureSetType mst = rat.getMeasureSet();
                Identifier variantIdIdentifier = new Identifier(IDENTIFIER_KEY_VARIATION, mst.getID().toString());
                List<Identifier> foundVariationIdentifiers = hearsayDAOBeanService.getIdentifierDAO().findByExample(variantIdIdentifier);
                if (CollectionUtils.isEmpty(foundVariationIdentifiers)) {
                    hearsayDAOBeanService.getIdentifierDAO().save(variantIdIdentifier);
                }

                List<Measure> measures = mst.getMeasure();
                if (CollectionUtils.isEmpty(measures)) {
                    return;
                }

                for (Measure measure : measures) {

                    if (!allowedTypes.contains(measure.getType())) {
                        continue;
                    }

                    List<XrefType> measureXrefs = measure.getXRef();

                    if (CollectionUtils.isNotEmpty(measureXrefs)) {
                        for (XrefType xref : measureXrefs) {
                            if ("dbSNP".equalsIgnoreCase(xref.getDB()) && "rs".equalsIgnoreCase(xref.getType())) {
                                String dbSNPId = xref.getID();
                                Identifier identifier = new Identifier(IDENTIFIER_KEY_SNP, String.format("rs%s", dbSNPId));
                                List<Identifier> foundSNPIdentifiers = hearsayDAOBeanService.getIdentifierDAO().findByExample(identifier);
                                if (CollectionUtils.isEmpty(foundSNPIdentifiers)) {
                                    hearsayDAOBeanService.getIdentifierDAO().save(identifier);
                                }
                            }
                        }
                    }

                    List<AttributeSet> attributeSetList = measure.getAttributeSet();

                    if (CollectionUtils.isNotEmpty(attributeSetList)) {

                        for (AttributeSet attributeSet : attributeSetList) {
                            List<XrefType> attributeXrefs = attributeSet.getXRef();

                            if (CollectionUtils.isNotEmpty(attributeXrefs)) {
                                for (XrefType xref : attributeXrefs) {
                                    if ("dbSNP".equalsIgnoreCase(xref.getDB()) && "rs".equalsIgnoreCase(xref.getType())) {
                                        String dbSNPId = xref.getID();
                                        Identifier identifier = new Identifier(IDENTIFIER_KEY_SNP, String.format("rs%s", dbSNPId));
                                        List<Identifier> foundSNPIdentifiers = hearsayDAOBeanService.getIdentifierDAO()
                                                .findByExample(identifier);
                                        if (CollectionUtils.isEmpty(foundSNPIdentifiers)) {
                                            hearsayDAOBeanService.getIdentifierDAO().save(identifier);
                                        }
                                    }
                                }
                            }
                        }

                    }

                }

            } catch (Exception e) {
                logger.error("Error", e);
                e.printStackTrace();
            }

        }

    }

    private void persistCanonicalAlleles(List<PublicSetType> publicSetTypeList) {
        try {
            ExecutorService es = Executors.newFixedThreadPool(4);
            for (PublicSetType pst : publicSetTypeList) {

                es.submit(() -> {

                    try {
                        ReferenceAssertionType rat = pst.getReferenceClinVarAssertion();
                        ReferenceAssertionType.ClinVarAccession clinVarAccession = rat.getClinVarAccession();
                        MeasureSetType mst = rat.getMeasureSet();

                        String prefix = pst.getTitle().substring(0, 3);
                        ReferenceSequenceType refSeqType = null;
                        for (ReferenceSequenceType referenceSequenceType : ReferenceSequenceType.values()) {
                            if (referenceSequenceType.getPrefixes().contains(prefix)) {
                                refSeqType = referenceSequenceType;
                                break;
                            }
                        }

                        CanonicalAlleleType canonicalAlleleType = null;
                        if (refSeqType != null) {
                            switch (refSeqType) {
                                case GENOMIC:
                                case RNA:
                                case TRANSCRIPT:
                                    canonicalAlleleType = CanonicalAlleleType.NUCLEOTIDE;
                                    break;
                                case PROTEIN:
                                    canonicalAlleleType = CanonicalAlleleType.AMINO_ACID;
                                    break;
                            }
                        }

                        if (canonicalAlleleType != null) {
                            CanonicalAllele canonicalAllele = new CanonicalAllele();
                            canonicalAllele.setActive("current".equals(rat.getRecordStatus()));
                            canonicalAllele.setVersion(clinVarAccession.getVersion().toString());
                            // TODO is this the right way to determine complexity???
                            canonicalAllele.setComplexityType(mst.getMeasure().size() > 1 ? ComplexityType.COMPLEX : ComplexityType.SIMPLE);
                            canonicalAllele.setType(canonicalAlleleType);

                            List<Identifier> foundIdentifiers = hearsayDAOBeanService.getIdentifierDAO()
                                    .findByExample(new Identifier(IDENTIFIER_KEY_VARIATION, mst.getID().toString()));
                            if (CollectionUtils.isNotEmpty(foundIdentifiers)) {
                                canonicalAllele.getIdentifiers().add(foundIdentifiers.get(0));
                            }

                            canonicalAllele.setId(hearsayDAOBeanService.getCanonicalAlleleDAO().save(canonicalAllele));
                        }

                    } catch (Exception e) {
                        logger.error("Error", e);
                        e.printStackTrace();
                    }
                });

            }
            es.shutdown();
            es.awaitTermination(20L, TimeUnit.MINUTES);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    private void persistMeasureAttributeContextualAlleles(List<PublicSetType> publicSetTypeList) {
        try {
            ExecutorService es = Executors.newFixedThreadPool(4);
            for (PublicSetType pst : publicSetTypeList) {

                es.submit(() -> {

                    try {
                        ReferenceAssertionType rat = pst.getReferenceClinVarAssertion();
                        MeasureSetType mst = rat.getMeasureSet();

                        List<CanonicalAllele> foundCanonicalAlleles = hearsayDAOBeanService.getCanonicalAlleleDAO()
                                .findByIdentifierSystemAndValue(IDENTIFIER_KEY_VARIATION, mst.getID().toString());
                        if (CollectionUtils.isEmpty(foundCanonicalAlleles)) {
                            logger.warn("Could not find CanonicalAllele: {}", mst.getID().toString());
                            return;
                        }
                        CanonicalAllele canonicalAllele = foundCanonicalAlleles.get(0);

                        List<Measure> measures = mst.getMeasure();
                        if (CollectionUtils.isNotEmpty(measures)) {

                            for (Measure measure : measures) {

                                if (!allowedTypes.contains(measure.getType())) {
                                    continue;
                                }

                                List<AttributeSet> attributeSetList = measure.getAttributeSet();

                                if (CollectionUtils.isEmpty(attributeSetList)) {
                                    continue;
                                }

                                List<XrefType> xrefs = measure.getXRef();

                                Identifier snpIdentifier = null;
                                if (CollectionUtils.isNotEmpty(xrefs)) {
                                    for (XrefType xref : xrefs) {
                                        if ("dbSNP".equalsIgnoreCase(xref.getDB()) && "rs".equalsIgnoreCase(xref.getType())) {
                                            List<Identifier> foundSNPIdentifiers = hearsayDAOBeanService.getIdentifierDAO()
                                                    .findByExample(new Identifier(IDENTIFIER_KEY_SNP, String.format("rs%s", xref.getID())));
                                            if (CollectionUtils.isNotEmpty(foundSNPIdentifiers)) {
                                                snpIdentifier = foundSNPIdentifiers.get(0);
                                            }
                                        }
                                    }
                                }

                                for (AttributeSet attributeSet : attributeSetList) {
                                    Attribute attribute = attributeSet.getAttribute();
                                    String attributeValue = attribute.getValue();
                                    String attributeType = attribute.getType();
                                    logger.info("attributeValue: {}, attributeType: {}", attributeValue, attributeType);
                                    ContextualAlleleNameType nameType = null;
                                    if (attributeValue.contains(":c.")) {
                                        nameType = ContextualAlleleNameType.HGVS_CDNA;
                                    } else if (attributeValue.contains(":g.")) {
                                        nameType = ContextualAlleleNameType.HGVS_GENOMIC;
                                    } else if (attributeValue.contains(":m.")) {
                                        nameType = ContextualAlleleNameType.HGVS_MITO;
                                    } else if (attributeValue.contains(":n.")) {
                                        nameType = ContextualAlleleNameType.HGVS_NCRNA;
                                    } else if (attributeValue.contains(":r.")) {
                                        nameType = ContextualAlleleNameType.HGVS_RNA;
                                    } else if (attributeValue.contains(":p.")) {
                                        nameType = ContextualAlleleNameType.HGVS_PROTEIN_1;
                                    } else {
                                        nameType = ContextualAlleleNameType.CUSTOM;
                                    }

                                    ReferenceSequence referenceSequence = null;
                                    switch (attributeType) {
                                        case "HGVS, genomic, top level":
                                            break;
                                        case "HGVS, coding, RefSeq":

                                            DNAVariantMutation variantMutation = HGVSParser.getInstance().parseDNAMutation(attributeValue);
                                            logger.info(variantMutation.toString());
                                            DNAChangeType changeType = variantMutation.getChangeType();
                                            if (changeType == null) {
                                                logger.warn("changeType is null: {}", attributeValue);
                                                return;
                                            }

                                            List<ReferenceSequence> foundReferenceSequences = hearsayDAOBeanService
                                                    .getReferenceSequenceDAO()
                                                    .findByIdentifierSystemAndValue(IDENTIFIER_KEY_NUCCORE, variantMutation.getAccession());

                                            if (CollectionUtils.isNotEmpty(foundReferenceSequences)) {
                                                referenceSequence = foundReferenceSequences.get(0);
                                            }

                                            if (referenceSequence == null) {
                                                logger.warn("No ReferenceSequences found: {}", attributeValue);
                                                return;
                                            }

                                            logger.debug("{}:{}", attributeType, attributeValue);
                                            AlleleInfo alleleInfo = variantMutation.getAlleleInfo();

                                            // TODO implement other AlleleInfo instances (Deletion, Insertion, etc.)
                                            if (alleleInfo instanceof SubstitutionAlleleInfo) {

                                                SubstitutionAlleleInfo substitutionAlleleInfo = (SubstitutionAlleleInfo) alleleInfo;

                                                ReferenceCoordinate referenceCoordinate = new ReferenceCoordinate();
                                                referenceCoordinate.setReferenceSequence(referenceSequence);
                                                referenceCoordinate.getIdentifiers().add(snpIdentifier);
                                                referenceCoordinate.setRefAllele(substitutionAlleleInfo.getWildtype());
                                                referenceCoordinate
                                                        .setId(hearsayDAOBeanService.getReferenceCoordinateDAO().save(referenceCoordinate));

                                                ContextualAllele contextualAllele = new ContextualAllele();
                                                contextualAllele.setCanonicalAllele(canonicalAllele);
                                                contextualAllele.setType(ContextualAlleleType.TRANSCRIPT);
                                                contextualAllele.setReferenceCoordinate(referenceCoordinate);
                                                contextualAllele.setAllele(substitutionAlleleInfo.getMutation());
                                                contextualAllele
                                                        .setId(hearsayDAOBeanService.getContextualAlleleDAO().save(contextualAllele));

                                                ContextualAlleleName contextualAlleleName = new ContextualAlleleName(attributeValue,
                                                        nameType);
                                                contextualAlleleName.setContextualAllele(contextualAllele);
                                                contextualAlleleName.setId(
                                                        hearsayDAOBeanService.getContextualAlleleNameDAO().save(contextualAlleleName));

                                                String location = substitutionAlleleInfo.getLocation();

                                                Pattern p = Pattern.compile("(\\d+)");
                                                Matcher m = p.matcher(location);
                                                if (m.matches()) {

                                                    Integer index = Integer.valueOf(m.group(1));

                                                    ExternalOffsetPosition endPosition = new ExternalOffsetPosition();
                                                    endPosition.setIndex(index);
                                                    endPosition
                                                            .setId(hearsayDAOBeanService.getExternalOffsetPositionDAO().save(endPosition));
                                                    referenceCoordinate.setEnd(endPosition);

                                                    ExternalOffsetPosition startPosition = new ExternalOffsetPosition();
                                                    startPosition.setIndex(index);
                                                    startPosition.setId(
                                                            hearsayDAOBeanService.getExternalOffsetPositionDAO().save(startPosition));
                                                    referenceCoordinate.setStart(startPosition);

                                                }

                                                p = Pattern.compile("(\\d+)([-|+])(\\d+)");
                                                m = p.matcher(location);
                                                if (m.matches()) {

                                                    DirectionType directionType = null;
                                                    for (DirectionType dt : DirectionType.values()) {
                                                        if (dt.getValue().equals(m.group(2))) {
                                                            directionType = dt;
                                                            break;
                                                        }
                                                    }
                                                    if (directionType != null) {

                                                        Integer index = Integer.valueOf(m.group(1));
                                                        Integer length = Integer.valueOf(m.group(3));

                                                        ExternalOffsetPosition startPosition = new ExternalOffsetPosition();
                                                        startPosition.setIndex(index);
                                                        startPosition.setDirection(directionType);
                                                        startPosition.setLength(length - 1);
                                                        startPosition.setId(
                                                                hearsayDAOBeanService.getExternalOffsetPositionDAO().save(startPosition));
                                                        referenceCoordinate.setStart(startPosition);

                                                        ExternalOffsetPosition endPosition = new ExternalOffsetPosition();
                                                        endPosition.setIndex(index);
                                                        endPosition.setDirection(directionType);
                                                        endPosition.setLength(length);
                                                        endPosition.setId(
                                                                hearsayDAOBeanService.getExternalOffsetPositionDAO().save(endPosition));
                                                        referenceCoordinate.setEnd(endPosition);

                                                    }

                                                }
                                                hearsayDAOBeanService.getReferenceCoordinateDAO().save(referenceCoordinate);

                                            }

                                            break;
                                        case "HGVS, RNA":
                                        case "HGVS, protein":
                                        case "HGVS, protein, RefSeq":
                                            break;
                                    }

                                }
                            }
                        }
                    } catch (Exception e) {
                        logger.error("Error", e);
                        e.printStackTrace();
                    }

                });
            }
            es.shutdown();
            es.awaitTermination(4L, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    private void persistMeasureNameContextualAlleles(List<PublicSetType> publicSetTypeList) {
        try {
            ExecutorService es = Executors.newFixedThreadPool(4);
            for (PublicSetType pst : publicSetTypeList) {

                es.submit(() -> {

                    try {
                        ReferenceAssertionType rat = pst.getReferenceClinVarAssertion();
                        MeasureSetType mst = rat.getMeasureSet();

                        List<CanonicalAllele> foundCanonicalAlleles = hearsayDAOBeanService.getCanonicalAlleleDAO()
                                .findByIdentifierSystemAndValue(IDENTIFIER_KEY_VARIATION, mst.getID().toString());
                        if (CollectionUtils.isEmpty(foundCanonicalAlleles)) {
                            logger.warn("Could not find CanonicalAllele: {}", mst.getID().toString());
                            return;
                        }
                        CanonicalAllele canonicalAllele = foundCanonicalAlleles.get(0);

                        List<Measure> measures = mst.getMeasure();
                        if (CollectionUtils.isNotEmpty(measures)) {

                            for (Measure measure : measures) {

                                if (!allowedTypes.contains(measure.getType())) {
                                    continue;
                                }

                                List<SetElementSetType> names = measure.getName();

                                if (CollectionUtils.isEmpty(names)) {
                                    continue;
                                }

                                List<ContextualAlleleName> contextualAlleleNameList = new ArrayList<ContextualAlleleName>();
                                for (SetElementSetType name : names) {
                                    ElementValue elementValue = name.getElementValue();
                                    String value = elementValue.getValue();
                                    ContextualAlleleNameType nameType = null;
                                    if (value.contains(":c.")) {
                                        nameType = ContextualAlleleNameType.HGVS_CDNA;
                                    } else if (value.contains(":g.")) {
                                        nameType = ContextualAlleleNameType.HGVS_GENOMIC;
                                    } else if (value.contains(":m.")) {
                                        nameType = ContextualAlleleNameType.HGVS_MITO;
                                    } else if (value.contains(":n.")) {
                                        nameType = ContextualAlleleNameType.HGVS_NCRNA;
                                    } else if (value.contains(":r.")) {
                                        nameType = ContextualAlleleNameType.HGVS_RNA;
                                    } else if (value.contains(":p.")) {
                                        nameType = ContextualAlleleNameType.HGVS_PROTEIN_1;
                                    } else {
                                        nameType = ContextualAlleleNameType.CUSTOM;
                                    }

                                    ContextualAlleleName contextualAlleleName = new ContextualAlleleName(elementValue.getValue(), nameType);
                                    contextualAlleleName
                                            .setPreferred("Preferred".equals(elementValue.getType()) ? Boolean.TRUE : Boolean.FALSE);
                                    contextualAlleleNameList.add(contextualAlleleName);
                                }

                                ContextualAllele contextualAllele = new ContextualAllele();
                                contextualAllele.setCanonicalAllele(canonicalAllele);
                                contextualAllele.setId(hearsayDAOBeanService.getContextualAlleleDAO().save(contextualAllele));

                                ContextualAlleleName preferredContextualAlleleName = null;
                                for (ContextualAlleleName contextualAlleleName : contextualAlleleNameList) {
                                    contextualAlleleName.setContextualAllele(contextualAllele);
                                    contextualAlleleName
                                            .setId(hearsayDAOBeanService.getContextualAlleleNameDAO().save(contextualAlleleName));
                                    if (contextualAlleleName.getPreferred()) {
                                        preferredContextualAlleleName = contextualAlleleName;
                                    }
                                }

                                if (preferredContextualAlleleName != null
                                        && !preferredContextualAlleleName.getType().equals(ContextualAlleleNameType.CUSTOM)) {
                                }

                            }
                        }
                    } catch (Exception e) {
                        logger.error("Error", e);
                        e.printStackTrace();
                    }

                });
            }
            es.shutdown();
            es.awaitTermination(4L, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

}
