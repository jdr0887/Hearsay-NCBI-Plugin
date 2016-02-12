package org.renci.hearsay.commands.ncbi;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Unmarshaller;

import org.apache.commons.collections.CollectionUtils;
import org.junit.Test;
import org.renci.clinvar.MeasureSetType;
import org.renci.clinvar.MeasureSetType.Measure;
import org.renci.clinvar.MeasureSetType.Measure.AttributeSet;
import org.renci.clinvar.MeasureSetType.Measure.AttributeSet.Attribute;
import org.renci.clinvar.PublicSetType;
import org.renci.clinvar.ReferenceAssertionType;
import org.renci.clinvar.ReleaseType;
import org.renci.clinvar.SetElementSetType;
import org.renci.hearsay.commands.ncbi.util.FTPUtil;

public class ParseClinVarTest {

    private static final List<String> allowedTypes = Arrays.asList("single nucleotide variant", "Duplication", "Deletion", "Indel",
            "inversion");

    @Test
    public void parseClinVar() {
        try {
            File clinvarDownload = FTPUtil.ncbiDownload("/pub/clinvar/xml", "ClinVarFullRelease_00-latest.xml.gz");
            JAXBContext jc = JAXBContext.newInstance(ReleaseType.class);
            Unmarshaller u = jc.createUnmarshaller();
            ReleaseType releaseType = (ReleaseType) u.unmarshal(new GZIPInputStream(new FileInputStream(clinvarDownload)));
            List<PublicSetType> publicSetType = releaseType.getClinVarSet();

            for (PublicSetType pst : publicSetType) {
                ReferenceAssertionType rat = pst.getReferenceClinVarAssertion();
                MeasureSetType mst = rat.getMeasureSet();

                if ("Variant".equals(mst.getType())) {

                    int count = 0;
                    for (Measure measure : mst.getMeasure()) {

                        if (!allowedTypes.contains(measure.getType())) {
                            continue;
                        }

                        List<AttributeSet> attributeSetList = measure.getAttributeSet();
                        for (AttributeSet attributeSet : attributeSetList) {
                            Attribute attribute = attributeSet.getAttribute();
                            String attributeValue = attribute.getValue();
                            String attributeType = attribute.getType();
                            if ("HGVS, genomic, top level".equals(attributeType)) {
                                count++;
                            }
                        }
                    }
                    if (count == 0) {
                        List<SetElementSetType> names = mst.getName();
                        if (CollectionUtils.isEmpty(names)) {
                            continue;
                        }
                        for (SetElementSetType name : names) {
                            System.out.println(name.getElementValue().getValue());
                        }
                    }
                    System.out.printf("'HGVS, genomic, top level' per Measures: %d/%d%n", count, mst.getMeasure().size());
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Test
    public void countVariantMeasuresPerMeasureSet() {
        try {
            File clinvarDownload = FTPUtil.ncbiDownload("/pub/clinvar/xml", "ClinVarFullRelease_00-latest.xml.gz");
            JAXBContext jc = JAXBContext.newInstance(ReleaseType.class);
            Unmarshaller u = jc.createUnmarshaller();
            ReleaseType releaseType = (ReleaseType) u.unmarshal(new GZIPInputStream(new FileInputStream(clinvarDownload)));
            List<PublicSetType> publicSetType = releaseType.getClinVarSet();
            publicSetType.forEach(a -> {
                ReferenceAssertionType rat = a.getReferenceClinVarAssertion();
                MeasureSetType mst = rat.getMeasureSet();
                System.out.println(mst.getMeasure().size());
            });
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Test
    public void printMeasureName() {
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(new File("/tmp", "types.txt")))) {
            File clinvarDownload = FTPUtil.ncbiDownload("/pub/clinvar/xml", "ClinVarFullRelease_00-latest.xml.gz");
            JAXBContext jc = JAXBContext.newInstance(ReleaseType.class);
            Unmarshaller u = jc.createUnmarshaller();
            ReleaseType releaseType = (ReleaseType) u.unmarshal(new GZIPInputStream(new FileInputStream(clinvarDownload)));
            List<PublicSetType> publicSetType = releaseType.getClinVarSet();

            List<String> allowedTypes = Arrays.asList("single nucleotide variant", "Duplication", "Deletion", "Indel", "inversion");

            // Set<String> referenceSequenceAccessionSet = new HashSet<String>();
            publicSetType.forEach(a -> {
                try {
                    ReferenceAssertionType rat = a.getReferenceClinVarAssertion();
                    MeasureSetType mst = rat.getMeasureSet();

                    if ("Variant".equals(mst.getType())) {
                        for (Measure measure : mst.getMeasure()) {
                            List<SetElementSetType> names = measure.getName();
                            if (CollectionUtils.isEmpty(names)) {
                                continue;
                            }
                            for (SetElementSetType name : names) {

                                if (!allowedTypes.contains(measure.getType())) {
                                    continue;
                                }

                                if ("Preferred".equals(name.getElementValue().getType())) {
                                    bw.write(String.format("%s:%s%n", measure.getType(), name.getElementValue().getValue()));
                                    bw.flush();
                                }
                            }
                            // String name = ev.getValue();
                            // int colonIndex = name.indexOf(":");
                            // if (colonIndex != -1) {
                            // name = name.substring(0, colonIndex);
                            // int startParenIndex = name.indexOf("(");
                            // if (startParenIndex != -1) {
                            // name = name.substring(0, startParenIndex);
                            // referenceSequenceAccessionSet.add(name);
                            // }
                            // }
                        }
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }
            });
            // referenceSequenceAccessionSet.forEach(a -> System.out.println(a));

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Test
    public void writeHGVS() {
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(new File("/tmp", "hgvs.txt")))) {
            File clinvarDownload = FTPUtil.ncbiDownload("/pub/clinvar/xml", "ClinVarFullRelease_00-latest.xml.gz");
            JAXBContext jc = JAXBContext.newInstance(ReleaseType.class);
            Unmarshaller u = jc.createUnmarshaller();
            ReleaseType releaseType = (ReleaseType) u.unmarshal(new GZIPInputStream(new FileInputStream(clinvarDownload)));
            List<PublicSetType> publicSetType = releaseType.getClinVarSet();

            // Set<String> referenceSequenceAccessionSet = new HashSet<String>();
            publicSetType.forEach(a -> {
                try {
                    ReferenceAssertionType rat = a.getReferenceClinVarAssertion();
                    MeasureSetType mst = rat.getMeasureSet();

                    if ("Variant".equals(mst.getType())) {
                        for (Measure measure : mst.getMeasure()) {

                            if (!allowedTypes.contains(measure.getType())) {
                                continue;
                            }

                            List<AttributeSet> attributeSetList = measure.getAttributeSet();
                            for (AttributeSet attributeSet : attributeSetList) {
                                Attribute attribute = attributeSet.getAttribute();
                                String attributeValue = attribute.getValue();
                                String attributeType = attribute.getType();

                                switch (attributeType) {
                                    case "HGVS, genomic, top level":
                                        break;
                                    case "HGVS, coding, RefSeq":
                                        bw.write(String.format("%s%n", attributeValue));
                                        bw.flush();
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
                    e.printStackTrace();
                }
            });

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Test
    public void printMolecularConsequenceValue() {
        try {
            File clinvarDownload = FTPUtil.ncbiDownload("/pub/clinvar/xml", "ClinVarFullRelease_00-latest.xml.gz");
            JAXBContext jc = JAXBContext.newInstance(ReleaseType.class);
            Unmarshaller u = jc.createUnmarshaller();
            ReleaseType releaseType = (ReleaseType) u.unmarshal(new GZIPInputStream(new FileInputStream(clinvarDownload)));
            List<PublicSetType> publicSetType = releaseType.getClinVarSet();
            Set<String> molecularConsequenceValueSet = new HashSet<String>();
            publicSetType.forEach(a -> {
                ReferenceAssertionType rat = a.getReferenceClinVarAssertion();
                MeasureSetType mst = rat.getMeasureSet();
                for (Measure measure : mst.getMeasure()) {

                    List<AttributeSet> attributeSetList = measure.getAttributeSet();
                    for (AttributeSet attributeSet : attributeSetList) {
                        Attribute attribute = attributeSet.getAttribute();
                        String attributeValue = attribute.getValue();
                        String attributeType = attribute.getType();

                        if ("MolecularConsequence".equals(attributeType)) {
                            molecularConsequenceValueSet.add(attributeValue);
                        }
                    }
                }
            });
            molecularConsequenceValueSet.forEach(a -> System.out.println(a));

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}
