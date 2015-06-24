package org.renci.hearsay.commands.ncbi;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.junit.Test;

public class Scratch {

    @Test
    public void test() {
        // String asdf = "join(1..3,4..5)";
        String asdf = "order(366..368,375..377,381..383,390..392,399..401)";

        Pattern p = Pattern.compile("^(join|order)\\((.+)\\)$");
        Matcher m = p.matcher(asdf);
        System.out.println(m.find());
        System.out.println(m.group(2));
        // System.out.println(m.group(2));
    }
}
