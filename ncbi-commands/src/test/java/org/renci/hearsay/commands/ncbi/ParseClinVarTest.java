package org.renci.hearsay.commands.ncbi;

import java.io.File;
import java.io.FileInputStream;
import java.util.List;
import java.util.zip.GZIPInputStream;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Unmarshaller;

import org.junit.Test;
import org.renci.clinvar.PublicSetType;
import org.renci.clinvar.ReferenceAssertionType;
import org.renci.clinvar.ReleaseType;
import org.renci.hearsay.commands.ncbi.util.FTPUtil;

public class ParseClinVarTest {

    @Test
    public void parseClinvRar() {
        File clinvarDownload = FTPUtil.ncbiDownload("/pub/clinvar/xml", "ClinVarFullRelease_00-latest.xml.gz");
        long start = System.currentTimeMillis();
        try {
            JAXBContext jc = JAXBContext.newInstance(ReleaseType.class);
            Unmarshaller u = jc.createUnmarshaller();
            ReleaseType releaseType = (ReleaseType) u.unmarshal(new GZIPInputStream(
                    new FileInputStream(clinvarDownload)));
            List<PublicSetType> publicSetType = releaseType.getClinVarSet();
            int count = 0;
            for (PublicSetType pst : publicSetType) {
                ReferenceAssertionType rat = pst.getReferenceClinVarAssertion();
                ReferenceAssertionType.ClinVarAccession clinVarAccession = rat.getClinVarAccession();
                count++;
            }
            System.out.println(count);
        } catch (Exception e) {
            e.printStackTrace();
        }
        long end = System.currentTimeMillis();
        System.out.printf("duration: %d seconds", (end - start) / 1000);

    }

}
