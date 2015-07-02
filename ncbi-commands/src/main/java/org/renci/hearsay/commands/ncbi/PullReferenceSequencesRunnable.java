package org.renci.hearsay.commands.ncbi;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang.math.IntRange;
import org.apache.commons.lang3.math.NumberUtils;
import org.renci.gbff.GBFFFilter;
import org.renci.gbff.GBFFManager;
import org.renci.gbff.filter.GBFFAndFilter;
import org.renci.gbff.filter.GBFFFeatureSourceOrganismNameFilter;
import org.renci.gbff.filter.GBFFFeatureTypeNameFilter;
import org.renci.gbff.filter.GBFFSequenceAccessionPrefixFilter;
import org.renci.gbff.filter.GBFFSourceOrganismNameFilter;
import org.renci.gbff.model.Feature;
import org.renci.gbff.model.Sequence;
import org.renci.gene2accession.G2AFilter;
import org.renci.gene2accession.G2AParser;
import org.renci.gene2accession.filter.G2AAndFilter;
import org.renci.gene2accession.filter.G2AAssemblyFilter;
import org.renci.gene2accession.filter.G2ATaxonIdFilter;
import org.renci.gene2accession.model.OrientationType;
import org.renci.gene2accession.model.Record;
import org.renci.hearsay.commands.ncbi.util.FTPUtil;
import org.renci.hearsay.dao.HearsayDAOBean;
import org.renci.hearsay.dao.HearsayDAOException;
import org.renci.hearsay.dao.model.Alignment;
import org.renci.hearsay.dao.model.Gene;
import org.renci.hearsay.dao.model.GenomeReference;
import org.renci.hearsay.dao.model.Identifier;
import org.renci.hearsay.dao.model.Location;
import org.renci.hearsay.dao.model.ReferenceSequence;
import org.renci.hearsay.dao.model.ReferenceSequenceChromosomeType;
import org.renci.hearsay.dao.model.ReferenceSequenceType;
import org.renci.hearsay.dao.model.Region;
import org.renci.hearsay.dao.model.RegionType;
import org.renci.hearsay.dao.model.StrandType;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class PullReferenceSequencesRunnable implements Runnable {

    private final Logger logger = LoggerFactory.getLogger(PullReferenceSequencesRunnable.class);

    private HearsayDAOBean hearsayDAOBean;

    public PullReferenceSequencesRunnable() {
        super();
    }

    @Override
    public void run() {
        logger.debug("ENTERING run()");

        String genomeReferenceVersion = null;

        File readmeCurrentReleaseFile = FTPUtil.ncbiDownload("/genomes/H_sapiens", "README_CURRENT_RELEASE");
        try (FileInputStream fis = new FileInputStream(readmeCurrentReleaseFile);
                InputStreamReader isr = new InputStreamReader(fis);
                BufferedReader br = new BufferedReader(isr)) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("ASSEMBLY NAME:")) {
                    genomeReferenceVersion = line.substring(14, line.length()).trim();
                    break;
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        logger.info("genomeReferenceVersion: {}", genomeReferenceVersion);
        // readmeCurrentReleaseFile.delete();

        // 5MB gzipped
        File refGeneFile = FTPUtil.ucscDownload("/goldenPath/hg38/database", "refGene.txt.gz");
        List<String> refGeneLines = new ArrayList<String>();
        try (FileInputStream fis = new FileInputStream(refGeneFile);
                GZIPInputStream gis = new GZIPInputStream(fis);
                InputStreamReader isr = new InputStreamReader(gis);
                BufferedReader br = new BufferedReader(isr)) {
            String line;
            while ((line = br.readLine()) != null) {
                refGeneLines.add(line);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        // refGeneFile.delete();

        // 380MB gzipped
        File genes2RefSeqFile = FTPUtil.ncbiDownload("/gene/DATA", "gene2refseq.gz");
        G2AParser gene2AccessionParser = G2AParser.getInstance();
        G2AAndFilter andFilter = new G2AAndFilter(Arrays.asList(new G2AFilter[] { new G2ATaxonIdFilter(9606),
                new G2AAssemblyFilter("Reference GRCh38.p2 Primary Assembly") }));
        List<Record> recordList = gene2AccessionParser.parse(andFilter, genes2RefSeqFile);

        Set<String> refseqVersionedAccessionSet = new HashSet<String>();
        for (Record record : recordList) {
            if (StringUtils.isNotEmpty(record.getRNANucleotideAccessionVersion())) {
                refseqVersionedAccessionSet.add(record.getRNANucleotideAccessionVersion());
            }
        }
        // genes2RefSeqFile.delete();

        List<String> accessionPrefixList = Arrays.asList(new String[] { "NM_", "NR_" });

        GBFFManager gbffMgr = GBFFManager.getInstance(1, true);

        List<GBFFFilter> filters = Arrays.asList(new GBFFFilter[] {
                new GBFFSequenceAccessionPrefixFilter(accessionPrefixList),
                new GBFFSourceOrganismNameFilter("Homo sapiens"),
                new GBFFFeatureSourceOrganismNameFilter("Homo sapiens"), new GBFFFeatureTypeNameFilter("CDS"),
                new GBFFFeatureTypeNameFilter("source") });

        List<File> vertebrateMammalianFileList = downloadVertebrateMammalianFiles();
        Pattern featureLocationPattern = Pattern.compile("^(join|order)\\((.+)\\)$");

        for (File f : vertebrateMammalianFileList) {

            try {
                List<Sequence> sequenceList = gbffMgr.deserialize(new GBFFAndFilter(filters), f);

                if (sequenceList != null && !sequenceList.isEmpty()) {

                    logger.info("sequenceList.size(): {}", sequenceList.size());

                    for (Sequence sequence : sequenceList) {

                        String refseqAccession = sequence.getAccession().contains(" ") ? sequence.getAccession()
                                .substring(0, sequence.getAccession().indexOf(" ")) : sequence.getAccession();
                        logger.info("refseqAccession: {}", refseqAccession);

                        String refSeqVersionedAccession = sequence.getVersion().trim().contains(" ") ? sequence
                                .getVersion().substring(0, sequence.getVersion().indexOf(" ")) : sequence.getVersion();
                        logger.info("refSeqVersionedAccession: {}", refSeqVersionedAccession);

                        if (!refseqVersionedAccessionSet.contains(refSeqVersionedAccession)) {
                            continue;
                        }

                        ReferenceSequence referenceSequence = new ReferenceSequence();

                        // set chromosome type
                        for (Feature feature : sequence.getFeatures()) {
                            if ("source".equals(feature.getType()) && feature.getQualifiers().containsKey("organism")) {
                                ReferenceSequenceChromosomeType chromosomeType = null;
                                for (ReferenceSequenceChromosomeType rsct : ReferenceSequenceChromosomeType.values()) {
                                    if (rsct.getValue().equals(feature.getQualifiers().getProperty("chromosome"))) {
                                        chromosomeType = rsct;
                                        break;
                                    }
                                }
                                referenceSequence.setChromosomeType(chromosomeType);
                                break;
                            }
                        }

                        // set type
                        for (Feature feature : sequence.getFeatures()) {
                            if ("CDS".equals(feature.getType())) {
                                referenceSequence.setType(ReferenceSequenceType.TRANSCRIPT);
                                break;
                            }
                        }

                        // set gene
                        for (Feature feature : sequence.getFeatures()) {
                            if ("gene".equals(feature.getType())) {
                                try {
                                    List<Gene> geneList = hearsayDAOBean.getGeneDAO().findBySymbol(
                                            feature.getQualifiers().getProperty("gene"));
                                    if (geneList != null && !geneList.isEmpty()) {
                                        referenceSequence.setGene(geneList.get(0));
                                    }
                                } catch (HearsayDAOException e) {
                                    e.printStackTrace();
                                }
                            }
                        }

                        // set GenomeReference
                        GenomeReference exampleGenomeReference = new GenomeReference();
                        exampleGenomeReference.setName(genomeReferenceVersion);
                        List<GenomeReference> potentialGenomeReferences = hearsayDAOBean.getGenomeReferenceDAO()
                                .findByExample(exampleGenomeReference);
                        if (potentialGenomeReferences != null && !potentialGenomeReferences.isEmpty()) {
                            referenceSequence.setGenomeReference(potentialGenomeReferences.get(0));
                        }

                        // initial save
                        referenceSequence.setId(hearsayDAOBean.getReferenceSequenceDAO().save(referenceSequence));
                        logger.info(referenceSequence.toString());

                        // set identifier
                        Identifier identifier = new Identifier("www.ncbi.nlm.nih.gov/nuccore", refseqAccession);
                        List<Identifier> possibleIdentifiers = hearsayDAOBean.getIdentifierDAO().findByExample(
                                identifier);
                        if (possibleIdentifiers != null && !possibleIdentifiers.isEmpty()) {
                            identifier = possibleIdentifiers.get(0);
                        } else {
                            identifier.setId(hearsayDAOBean.getIdentifierDAO().save(identifier));
                        }
                        logger.info(identifier.toString());
                        referenceSequence.getIdentifiers().add(identifier);

                        hearsayDAOBean.getReferenceSequenceDAO().save(referenceSequence);

                        // add alignments
                        final Alignment alignment = new Alignment();
                        alignment.setReferenceSequence(referenceSequence);
                        for (Record record : recordList) {
                            if (refSeqVersionedAccession.equals(record.getRNANucleotideAccessionVersion())) {
                                Location genomicLocation = new Location(record.getGenomicStartPosition(),
                                        record.getGenomicEndPosition());
                                genomicLocation.setId(hearsayDAOBean.getLocationDAO().save(genomicLocation));
                                alignment.setGenomicLocation(genomicLocation);
                                alignment
                                        .setStrandType(record.getOrientation().equals(OrientationType.MINUS) ? StrandType.MINUS
                                                : StrandType.PLUS);
                                break;
                            }
                        }

                        // add protein info to alignment
                        for (Feature feature : sequence.getFeatures()) {
                            if ("CDS".equals(feature.getType())) {
                                // String ccdsId = feature.getDBXrefs().getProperty("CCDS");
                                // if (StringUtils.isNotEmpty(ccdsId)) {
                                // alignment.setStrandType(ccdsExons2StrandMap.get(ccdsId));
                                // }
                                String range = feature.getLocation();
                                String start = range.substring(0, range.indexOf(".."));
                                String stop = range.substring(range.indexOf("..") + 2, range.length());
                                if (NumberUtils.isNumber(start) && NumberUtils.isNumber(stop)) {
                                    Location proteinLocation = new Location(Integer.valueOf(start),
                                            Integer.valueOf(stop));
                                    proteinLocation.setId(hearsayDAOBean.getLocationDAO().save(proteinLocation));
                                    alignment.setProteinLocation(proteinLocation);
                                }
                                alignment.setProtein(feature.getQualifiers().getProperty("protein_id"));
                                break;
                            }
                        }

                        alignment.setId(hearsayDAOBean.getAlignmentDAO().save(alignment));
                        logger.info(alignment.toString());

                        referenceSequence.getAlignments().add(alignment);
                        hearsayDAOBean.getReferenceSequenceDAO().save(referenceSequence);

                        // add exons to alignment
                        List<Region> regions = new ArrayList<Region>();
                        for (Feature feature : sequence.getFeatures()) {
                            if ("exon".equals(feature.getType())) {
                                Region region = new Region(RegionType.EXON);
                                String range = feature.getLocation();
                                String start = null;
                                String stop = null;
                                if (alignment.getStrandType().equals(StrandType.MINUS)) {
                                    stop = range.substring(0, range.indexOf(".."));
                                    start = range.substring(range.indexOf("..") + 2, range.length());
                                } else {
                                    start = range.substring(0, range.indexOf(".."));
                                    stop = range.substring(range.indexOf("..") + 2, range.length());
                                }

                                if (NumberUtils.isNumber(start) && NumberUtils.isNumber(stop)) {
                                    Location transcriptLocation = new Location(Integer.valueOf(start),
                                            Integer.valueOf(stop));
                                    transcriptLocation.setId(hearsayDAOBean.getLocationDAO().save(transcriptLocation));
                                    logger.debug(transcriptLocation.toString());
                                    region.setTranscriptLocation(transcriptLocation);
                                }
                                regions.add(region);
                            }
                        }

                        Collections.sort(regions, new Comparator<Region>() {
                            @Override
                            public int compare(Region r1, Region r2) {
                                if (alignment.getStrandType().equals(StrandType.MINUS)) {
                                    return r2.getTranscriptLocation().getStart()
                                            .compareTo(r1.getTranscriptLocation().getStart());
                                } else {
                                    return r1.getTranscriptLocation().getStart()
                                            .compareTo(r2.getTranscriptLocation().getStart());
                                }
                            }
                        });

                        // refseq accession = NM_182701
                        // genomic start = 28503295
                        // genomic end = 28515792
                        // exon starts = 28503295,28505702,28506311,28510750,28515656
                        // exon ends = 28504498,28505802,28506429,28510904,28515793

                        List<Location> genomicLocationList = new ArrayList<Location>();
                        for (String line : refGeneLines) {
                            Scanner scanner = new Scanner(line).useDelimiter("\t");
                            String bin = scanner.next();
                            String name = scanner.next();
                            if (!name.equals(refseqAccession)) {
                                continue;
                            }
                            String chromosome = scanner.next();
                            if (!chromosome.equals(String.format("chr%s", referenceSequence.getChromosomeType()
                                    .getValue()))) {
                                continue;
                            }
                            String strand = scanner.next();
                            if (!strand.equals(alignment.getStrandType().getValue())) {
                                continue;
                            }
                            String txStart = scanner.next();
                            String txEnd = scanner.next();
                            String cdsStart = scanner.next();
                            String cdsEnd = scanner.next();
                            String exonCount = scanner.next();

                            String exonStarts = scanner.next();
                            String[] exonStartsSplit = exonStarts.split(",");

                            String exonEnds = scanner.next();
                            String[] exonEndsSplit = exonEnds.split(",");

                            for (int i = 0; i < Integer.valueOf(exonCount); i++) {
                                Location genomicLocation = new Location(Integer.valueOf(exonStartsSplit[i]),
                                        Integer.valueOf(exonEndsSplit[i]));
                                logger.debug(genomicLocation.toString());
                                genomicLocationList.add(genomicLocation);
                            }
                        }

                        if (regions.size() != genomicLocationList.size()) {
                            logger.debug("regions.size(): {}", regions.size());
                            logger.debug("genomicLocationList.size(): {}", genomicLocationList.size());
                            logger.error("region size != genomic location size");
                            continue;
                        }

                        Collections.sort(genomicLocationList, new Comparator<Location>() {
                            @Override
                            public int compare(Location l1, Location l2) {
                                return l1.getStart().compareTo(l2.getStart());
                            }
                        });

                        // initial exon regions
                        for (int i = 0; i < regions.size(); i++) {
                            Region region = regions.get(i);
                            region.setAlignment(alignment);
                            Location regionLocation = genomicLocationList.get(i);
                            regionLocation.setId(hearsayDAOBean.getLocationDAO().save(regionLocation));
                            logger.debug(regionLocation.toString());
                            region.setRegionLocation(regionLocation);
                            region.setId(hearsayDAOBean.getRegionDAO().save(region));
                            logger.debug(region.toString());
                            alignment.getRegions().add(region);
                        }

                        hearsayDAOBean.getAlignmentDAO().save(alignment);

                        // adding utr regions
                        List<Region> utrRegionList = new ArrayList<Region>();
                        Iterator<Region> regionIter = alignment.getRegions().iterator();
                        while (regionIter.hasNext()) {
                            Region region = regionIter.next();

                            Region firstRegion = alignment.getRegions().get(0);
                            Region lastRegion = alignment.getRegions().get(alignment.getRegions().size() - 1);

                            if (alignment.getProteinLocation().getStart()
                                    .equals(firstRegion.getTranscriptLocation().getStart())
                                    && alignment.getProteinLocation().getStop()
                                            .equals(lastRegion.getTranscriptLocation().getStop())) {
                                continue;
                            }

                            Location regionLocation = region.getRegionLocation();
                            int regionStart = regionLocation.getStart();
                            int regionStop = regionLocation.getStop();

                            Location transcriptLocation = region.getTranscriptLocation();
                            IntRange transcriptRange = transcriptLocation.toRange();

                            int transcriptStart = transcriptLocation.getStart();
                            int transcriptStop = transcriptLocation.getStop();

                            if (alignment.getStrandType().equals(StrandType.MINUS)
                                    && transcriptRange.containsInteger(alignment.getProteinLocation().getStart())) {

                                transcriptLocation.setStop(alignment.getProteinLocation().getStart());
                                hearsayDAOBean.getLocationDAO().save(transcriptLocation);

                                int diff = transcriptStart - transcriptLocation.getStop();
                                regionLocation.setStop(regionLocation.getStart() + diff);
                                hearsayDAOBean.getLocationDAO().save(regionLocation);

                                Region newRegion = new Region(RegionType.UTR5);
                                newRegion.setAlignment(alignment);

                                Location newTranscriptLocation = new Location(
                                        alignment.getProteinLocation().getStart() - 1, transcriptStop);
                                newTranscriptLocation
                                        .setId(hearsayDAOBean.getLocationDAO().save(newTranscriptLocation));
                                newRegion.setTranscriptLocation(newTranscriptLocation);

                                Location newRegionLocation = new Location(regionLocation.getStop() + 1,
                                        regionLocation.getStop() + 1 + newRegion.getTranscriptLocation().diff());
                                newRegionLocation.setId(hearsayDAOBean.getLocationDAO().save(newRegionLocation));
                                newRegion.setRegionLocation(newRegionLocation);

                                newRegion.setId(hearsayDAOBean.getRegionDAO().save(newRegion));
                                utrRegionList.add(newRegion);
                            }

                            if (alignment.getStrandType().equals(StrandType.MINUS)
                                    && transcriptRange.containsInteger(alignment.getProteinLocation().getStop())) {

                                transcriptLocation.setStart(alignment.getProteinLocation().getStop());
                                hearsayDAOBean.getLocationDAO().save(transcriptLocation);

                                int diff = alignment.getProteinLocation().getStop() - transcriptStop;
                                regionLocation.setStart(regionLocation.getStop() - diff);
                                hearsayDAOBean.getLocationDAO().save(regionLocation);

                                Region newRegion = new Region(RegionType.UTR3);
                                newRegion.setAlignment(alignment);

                                Location newTranscriptLocation = new Location(transcriptStart, alignment
                                        .getProteinLocation().getStop() + 1);
                                newTranscriptLocation
                                        .setId(hearsayDAOBean.getLocationDAO().save(newTranscriptLocation));
                                newRegion.setTranscriptLocation(newTranscriptLocation);

                                Location newRegionLocation = new Location(regionLocation.getStart() - 2
                                        - newTranscriptLocation.diff(), regionLocation.getStart() - 1);
                                newRegionLocation.setId(hearsayDAOBean.getLocationDAO().save(newRegionLocation));
                                newRegion.setRegionLocation(newRegionLocation);

                                newRegion.setId(hearsayDAOBean.getRegionDAO().save(newRegion));
                                utrRegionList.add(newRegion);
                            }

                            if (alignment.getStrandType().equals(StrandType.PLUS)
                                    && transcriptRange.containsInteger(alignment.getProteinLocation().getStart())) {

                                transcriptLocation.setStart(alignment.getProteinLocation().getStart());
                                hearsayDAOBean.getLocationDAO().save(transcriptLocation);

                                regionLocation.setStart(regionStop - transcriptLocation.diff());
                                hearsayDAOBean.getLocationDAO().save(regionLocation);

                                Region newRegion = new Region(RegionType.UTR5);
                                newRegion.setAlignment(alignment);

                                Location newTranscriptLocation = new Location(transcriptStart, alignment
                                        .getProteinLocation().getStart() - 1);
                                newTranscriptLocation
                                        .setId(hearsayDAOBean.getLocationDAO().save(newTranscriptLocation));
                                newRegion.setTranscriptLocation(newTranscriptLocation);

                                Location newRegionLocation = new Location(regionLocation.getStart()
                                        - newTranscriptLocation.diff(), regionLocation.getStart() - 1);
                                newRegionLocation.setId(hearsayDAOBean.getLocationDAO().save(newRegionLocation));
                                newRegion.setRegionLocation(newRegionLocation);

                                newRegion.setId(hearsayDAOBean.getRegionDAO().save(newRegion));
                                utrRegionList.add(newRegion);

                            }

                            if (alignment.getStrandType().equals(StrandType.PLUS)
                                    && transcriptRange.containsInteger(alignment.getProteinLocation().getStop())) {

                                transcriptLocation.setStop(alignment.getProteinLocation().getStop());
                                hearsayDAOBean.getLocationDAO().save(transcriptLocation);

                                regionLocation.setStop(regionStart + transcriptLocation.diff());
                                hearsayDAOBean.getLocationDAO().save(regionLocation);

                                Region newRegion = new Region(RegionType.UTR3);
                                newRegion.setAlignment(alignment);

                                Location newTranscriptLocation = new Location(
                                        alignment.getProteinLocation().getStop() + 1, transcriptStop);
                                newTranscriptLocation
                                        .setId(hearsayDAOBean.getLocationDAO().save(newTranscriptLocation));
                                newRegion.setTranscriptLocation(newTranscriptLocation);

                                Location newRegionLocation = new Location(regionStop - newTranscriptLocation.diff(),
                                        regionStop);
                                newRegionLocation.setId(hearsayDAOBean.getLocationDAO().save(newRegionLocation));
                                newRegion.setRegionLocation(newRegionLocation);

                                newRegion.setId(hearsayDAOBean.getRegionDAO().save(newRegion));
                                utrRegionList.add(newRegion);

                            }

                        }

                        alignment.getRegions().addAll(utrRegionList);
                        hearsayDAOBean.getAlignmentDAO().save(alignment);

                        for (Region region : alignment.getRegions()) {

                            if (alignment.getStrandType().equals(StrandType.PLUS)
                                    && region.getTranscriptLocation().getStop() < alignment.getProteinLocation()
                                            .getStart()) {
                                region.setRegionType(RegionType.UTR5);
                            }

                            if (alignment.getStrandType().equals(StrandType.PLUS)
                                    && region.getTranscriptLocation().getStop() > alignment.getProteinLocation()
                                            .getStop()) {
                                region.setRegionType(RegionType.UTR3);
                            }

                            if (alignment.getStrandType().equals(StrandType.MINUS)
                                    && region.getTranscriptLocation().getStop() < alignment.getProteinLocation()
                                            .getStart()) {
                                region.setRegionType(RegionType.UTR5);
                            }

                            if (alignment.getStrandType().equals(StrandType.MINUS)
                                    && region.getTranscriptLocation().getStop() > alignment.getProteinLocation()
                                            .getStop()) {
                                region.setRegionType(RegionType.UTR3);
                            }

                            hearsayDAOBean.getRegionDAO().save(region);
                        }

                        // adding intron regions
                        List<Region> intronRegionList = new ArrayList<Region>();
                        regionIter = regions.iterator();
                        Region previousRegion = null;
                        while (regionIter.hasNext()) {
                            Region currentRegion = regionIter.next();
                            if (previousRegion == null) {
                                previousRegion = currentRegion;
                                continue;
                            }

                            if (!currentRegion.getRegionLocation().getStart()
                                    .equals(previousRegion.getRegionLocation().getStop() - 1)) {
                                Region region = new Region(RegionType.INTRON);
                                region.setAlignment(alignment);
                                Location regionLocation = new Location(
                                        previousRegion.getRegionLocation().getStop() + 1, currentRegion
                                                .getRegionLocation().getStart() - 1);
                                regionLocation.setId(hearsayDAOBean.getLocationDAO().save(regionLocation));
                                region.setRegionLocation(regionLocation);
                                region.setId(hearsayDAOBean.getRegionDAO().save(region));
                                intronRegionList.add(region);
                            }

                            previousRegion = currentRegion;
                        }

                        alignment.getRegions().addAll(intronRegionList);
                        hearsayDAOBean.getAlignmentDAO().save(alignment);

                        // add features
                        List<String> inclusionPatterns = Arrays.asList(new String[] { "misc_feature", "polyA_signal",
                                "polyA_site", "transit_peptide", "mat_peptide", "sig_peptide", "unsure", "stem_loop",
                                "protein_bind", "repeat_region", "prim_transcript", "proprotein", "LTR", "TATA_signal",
                                "primer_bind", "terminator", "misc_difference", "misc_binding", "RBS", "misc_signal",
                                "J_segment", "C_region", "conflict", "promoter", "ncRNA", "modified_base" });

                        for (Feature feature : sequence.getFeatures()) {

                            if (!inclusionPatterns.contains(feature.getType())) {
                                continue;
                            }
                            logger.debug(feature.toString());
                            String location = feature.getLocation();
                            org.renci.hearsay.dao.model.Feature hearsayFeature = new org.renci.hearsay.dao.model.Feature();
                            hearsayFeature.setReferenceSequence(referenceSequence);
                            hearsayFeature.setType(feature.getType());
                            String note = feature.getQualifiers().getProperty("note");

                            if (StringUtils.isNotEmpty(note)) {
                                hearsayFeature.setNote(note);
                            }

                            hearsayFeature.setId(hearsayDAOBean.getFeatureDAO().save(hearsayFeature));

                            if (location.startsWith("join") || location.startsWith("order")) {

                                Matcher m = featureLocationPattern.matcher(location);
                                m.find();

                                Scanner scanner = new Scanner(m.group(2)).useDelimiter(",");
                                while (scanner.hasNext()) {
                                    String range = scanner.next();
                                    String startValue = range.substring(0, range.indexOf(".."));
                                    String stopValue = range.substring(range.indexOf("..") + 2, range.length());
                                    if (NumberUtils.isNumber(startValue) && NumberUtils.isNumber(stopValue)) {
                                        Location l = new Location(Integer.valueOf(startValue),
                                                Integer.valueOf(stopValue));
                                        l.setId(hearsayDAOBean.getLocationDAO().save(l));
                                        hearsayFeature.getLocations().add(l);
                                    }

                                }

                            } else if (NumberUtils.isNumber(location)) {
                                Location l = new Location(Integer.valueOf(location), Integer.valueOf(location));
                                l.setId(hearsayDAOBean.getLocationDAO().save(l));
                                hearsayFeature.getLocations().add(l);
                            }

                            hearsayDAOBean.getFeatureDAO().save(hearsayFeature);
                            logger.info(hearsayFeature.toString());
                            referenceSequence.getFeatures().add(hearsayFeature);

                        }

                        hearsayDAOBean.getReferenceSequenceDAO().save(referenceSequence);
                        logger.info(referenceSequence.toString());

                    }

                }

                // f.delete();
            } catch (Exception e) {
                logger.error(e.getMessage(), e);
            }
        }

        logger.debug("FINISHED run()");
    }

    public List<File> downloadVertebrateMammalianFiles() {

        List<File> ret = new ArrayList<File>();

        File tmpDir = new File("/tmp");
        for (File f : tmpDir.listFiles()) {
            if (f.getName().startsWith("vertebrate_mammalian")) {
                ret.add(f);
            }
        }

        // FTPClient ftpClient = new FTPClient();
        // try {
        // ftpClient.connect("ftp.ncbi.nlm.nih.gov");
        //
        // ftpClient.login("anonymous", "anonymous");
        // ftpClient.setFileType(FTP.BINARY_FILE_TYPE);
        // ftpClient.enterLocalPassiveMode();
        //
        // int reply = ftpClient.getReplyCode();
        // if (!FTPReply.isPositiveCompletion(reply)) {
        // ftpClient.disconnect();
        // logger.error("FTP server refused connection.");
        // return null;
        // }
        //
        // List<FTPFile> ftpFileList = Arrays.asList(ftpClient.listFiles("/refseq/release/vertebrate_mammalian/",
        // new FTPFileFilter() {
        //
        // @Override
        // public boolean accept(FTPFile ftpFile) {
        // if (ftpFile != null && ftpFile.getName().endsWith("rna.gbff.gz")) {
        // return true;
        // }
        // return false;
        // }
        // }));
        //
        // for (FTPFile ftpFile : ftpFileList) {
        // File tmpFile = new File(System.getProperty("java.io.tmpdir", "/tmp"), ftpFile.getName());
        // try (OutputStream fos = new BufferedOutputStream(new FileOutputStream(tmpFile))) {
        // ftpClient.retrieveFile(String.format("/refseq/release/vertebrate_mammalian/%s", ftpFile.getName()),
        // fos);
        // fos.flush();
        // }
        // ret.add(tmpFile);
        // }
        //
        // } catch (IOException e) {
        // e.printStackTrace();
        // } finally {
        // try {
        // if (ftpClient.isConnected()) {
        // ftpClient.disconnect();
        // }
        // } catch (IOException e) {
        // e.printStackTrace();
        // }
        // }
        return ret;
    }

    public HearsayDAOBean getHearsayDAOBean() {
        return hearsayDAOBean;
    }

    public void setHearsayDAOBean(HearsayDAOBean hearsayDAOBean) {
        this.hearsayDAOBean = hearsayDAOBean;
    }

}
