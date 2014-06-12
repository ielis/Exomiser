/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package de.charite.compbio.exomiser.util;

import de.charite.compbio.exomiser.priority.PriorityType;
import de.charite.compbio.exomiser.util.ExomiserSettings.Builder;
import jannovar.common.ModeOfInheritance;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import static org.hamcrest.CoreMatchers.equalTo;
import static org.hamcrest.CoreMatchers.instanceOf;
import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.CoreMatchers.nullValue;
import static org.junit.Assert.*;
import org.junit.Before;
import org.junit.Test;

/**
 * Tests for {@code de.charite.compbio.exomiser.util.ExomiserSettings}.
 *
 * @author Jules Jacobsen <jules.jacobsen@sanger.ac.uk>
 */
public class ExomiserSettingsTest {

    Builder builder;

    //
    private static final Path VCF_PATH_NOT_SET = null;
    private static final Path VCF_PATH = Paths.get("test.vcf");
    private static final Path PED_PATH_NOT_SET = null;
    private static final Path PED_PATH = Paths.get("test.ped");
    private static final float MAXIMUM_FREQUENCY_DEFAULT = 100.0f;
    private static final float MAXIMUM_FREQUENCY = 42.24f;
    private static final float MIMIMUM_QUALITY_DEFAULT = 0.0f;
    private static final float MIMIMUM_QUALITY = 666.24f;
    private static final String GENETIC_INTERVAL_DEFAULT = "";
    private static final String GENETIC_INTERVAL = "chr2:12345-67890";
    private static final boolean INCLUDE_PATHOGENIC_DEFAULT = false;
    private static final boolean INCLUDE_PATHOGENIC = true;
    private static final boolean REMOVE_DBSNP_DEFAULT = false;
    private static final boolean REMOVE_DBSNP = true;
    private static final boolean REMOVE_OFF_TARGET_VARIANTS_DEFAULT = true;
    private static final boolean REMOVE_OFF_TARGET_VARIANTS = false;
    private static final String CANDIDATE_GENE_NAME_DEFAULT = "";
    private static final String CANDIDATE_GENE_NAME = "ADH1";
    private static final ModeOfInheritance MODE_OF_INHERITANCE = ModeOfInheritance.AUTOSOMAL_DOMINANT;
    private static final ModeOfInheritance MODE_OF_INHERITANCE_DEFAULT = ModeOfInheritance.UNINITIALIZED;
    private static final String DISEASE_STRING_DEFAULT = "";
    private static final String DISEASE_STRING = "OMIM:100100";
    private static final List<String> HPO_LIST_DEFAULT = new ArrayList<>();
    private static final List<String> HPO_LIST = new ArrayList<>(Arrays.asList("HPO:123456"));
    private static final List<Integer> SEED_GENE_LIST_DEFAULT = new ArrayList<>();
    private static final List<Integer> SEED_GENE_LIST = new ArrayList<>(Arrays.asList(1, 23, 56));
    private static final int NUMBER_OF_GENES_TO_SHOW_DEFAULT = 0;
    private static final int NUMBER_OF_GENES_TO_SHOW = 12438803;
    private static final String OUT_FILE_NAME_DEFAULT = "exomiser.html";
    private static final String OUT_FILE_NAME = "wibbler.tab";
    private static final OutputFormat OUTPUT_FORMAT_DEFAULT = OutputFormat.HTML;
    private static final OutputFormat OUTPUT_FORMAT = OutputFormat.TAB;

    public ExomiserSettingsTest() {
    }

    @Before
    public void setUp() {
        builder = new ExomiserSettings.Builder();
    }

    @Test
    public void testThatTheBuilderProducesDefaultExomiserSettingsObject() {
        ExomiserSettings settings = builder.build();
        assertThat(settings, instanceOf(ExomiserSettings.class));
        System.out.println(settings);
        assertThat(settings.getVcfPath(), equalTo(VCF_PATH_NOT_SET));
        assertThat(settings.getPedPath(), equalTo(PED_PATH_NOT_SET));
        assertThat(settings.getPrioritiserType(), equalTo(PriorityType.NOT_SET));
        assertThat(settings.getMaximumFrequency(), equalTo(MAXIMUM_FREQUENCY_DEFAULT));
        assertThat(settings.getMinimumQuality(), equalTo(MIMIMUM_QUALITY_DEFAULT));
        assertThat(settings.getGeneticInterval(), equalTo(GENETIC_INTERVAL_DEFAULT));
        assertThat(settings.includePathogenic(), is(INCLUDE_PATHOGENIC_DEFAULT));
        assertThat(settings.removeDbSnp(), is(REMOVE_DBSNP_DEFAULT));
        assertThat(settings.removeOffTargetVariants(), is(REMOVE_OFF_TARGET_VARIANTS_DEFAULT));
        assertThat(settings.getCandidateGene(), equalTo(CANDIDATE_GENE_NAME_DEFAULT));
        assertThat(settings.getModeOfInheritance(), equalTo(MODE_OF_INHERITANCE_DEFAULT));
        assertThat(settings.getDiseaseId(), equalTo(DISEASE_STRING_DEFAULT));
        assertThat(settings.getHpoIds(), equalTo(HPO_LIST_DEFAULT));
        assertThat(settings.getSeedGeneList(), equalTo(SEED_GENE_LIST_DEFAULT));
        assertThat(settings.getNumberOfGenesToShow(), equalTo(NUMBER_OF_GENES_TO_SHOW_DEFAULT));
        assertThat(settings.getOutFileName(), equalTo(OUT_FILE_NAME_DEFAULT));
        assertThat(settings.getOutputFormat(), equalTo(OUTPUT_FORMAT_DEFAULT));

    }

    /**
     * Test of getVcfPath method, of class ExomiserSettings.
     */
    @Test
    public void testThatGetVcfPathReturnsAPath() {

        builder.vcfFilePath(VCF_PATH);
        ExomiserSettings settings = builder.build();

        assertThat(settings.getVcfPath(), equalTo(VCF_PATH));
    }

    /**
     * Test of getVcfPath method, of class ExomiserSettings.
     */
    @Test
    public void testThatTheDefaultVcfPathIsNull() {
        ExomiserSettings settings = builder.build();
        assertThat(settings.getVcfPath(), nullValue());
    }

    /**
     * Test of getPedPath method, of class ExomiserSettings.
     */
    @Test
    public void testThatGetPedPathReturnsAPath() {

        builder.pedFilePath(PED_PATH);
        ExomiserSettings settings = builder.build();

        assertThat(settings.getPedPath(), equalTo(PED_PATH));
    }

    /**
     * Test of getPrioritiserType method, of class ExomiserSettings.
     */
    @Test
    public void testThatBuilderProducesAnUndefinedPriorityTypeAsDefault() {
        ExomiserSettings settings = builder.build();
        assertThat(settings.getPrioritiserType(), equalTo(PriorityType.NOT_SET));
    }

    /**
     * Test of getPrioritiserType method, of class ExomiserSettings.
     */
    @Test
    public void testThatBuilderProducesTheSpecifiedPriorityType() {
        builder.usePrioritiser(PriorityType.OMIM_PRIORITY);
        ExomiserSettings settings = builder.build();
        assertThat(settings.getPrioritiserType(), equalTo(PriorityType.OMIM_PRIORITY));
    }

    /**
     * Test of getMaximumFrequency method, of class ExomiserSettings.
     */
    @Test
    public void testThatBuilderProducesMaximumFrequencyDefault() {
        ExomiserSettings settings = builder.build();
        assertThat(settings.getMaximumFrequency(), equalTo(MAXIMUM_FREQUENCY_DEFAULT));
    }

    /**
     * Test of getMaximumFrequency method, of class ExomiserSettings.
     */
    @Test
    public void testThatBuilderProducesMaximumFrequencySpecified() {
        builder.maximumFrequency(MAXIMUM_FREQUENCY);
        ExomiserSettings settings = builder.build();
        assertThat(settings.getMaximumFrequency(), equalTo(MAXIMUM_FREQUENCY));
    }

    /**
     * Test of getMinimumQuality method, of class ExomiserSettings.
     */
    @Test
    public void testThatBuilderProducesDefaultMinimumQuality() {
        ExomiserSettings settings = builder.build();
        assertThat(settings.getMinimumQuality(), equalTo(MIMIMUM_QUALITY_DEFAULT));
    }

    @Test
    public void testThatBuilderProducesMinimumQualitySpecified() {
        builder.minimumQuality(MIMIMUM_QUALITY);
        ExomiserSettings settings = builder.build();
        assertThat(settings.getMinimumQuality(), equalTo(MIMIMUM_QUALITY));
    }

    /**
     * Test of getGeneticInterval method, of class ExomiserSettings.
     */
    @Test
    public void testThatBuilderProducesGeneticIntervalDefault() {
        ExomiserSettings settings = builder.build();
        assertThat(settings.getGeneticInterval(), equalTo(GENETIC_INTERVAL_DEFAULT));
    }

    @Test
    public void testThatBuilderProducesGeneticIntervalSpecified() {
        builder.geneticInterval(GENETIC_INTERVAL);
        ExomiserSettings settings = builder.build();
        assertThat(settings.getGeneticInterval(), equalTo(GENETIC_INTERVAL));
    }

    /**
     * Test of includePathogenic method, of class ExomiserSettings.
     */
    @Test
    public void testThatBuilderProducesIncludePathogenicDefault() {
        ExomiserSettings settings = builder.build();
        assertThat(settings.includePathogenic(), is(INCLUDE_PATHOGENIC_DEFAULT));
    }

    @Test
    public void testThatBuilderProducesIncludePathogenicWhenSet() {
        builder.includePathogenic(INCLUDE_PATHOGENIC);
        ExomiserSettings settings = builder.build();
        assertThat(settings.includePathogenic(), is(INCLUDE_PATHOGENIC));
    }

    /**
     * Test of removeDbSnp method, of class ExomiserSettings.
     */
    @Test
    public void testThatBuilderProducesRemoveDbSnpDefault() {
        ExomiserSettings settings = builder.build();
        assertThat(settings.removeDbSnp(), is(REMOVE_DBSNP_DEFAULT));
    }

    @Test
    public void testThatBuilderProducesRemoveDbSnpWhenSet() {
        builder.removeDbSnp(REMOVE_DBSNP);
        ExomiserSettings settings = builder.build();
        assertThat(settings.removeDbSnp(), is(REMOVE_DBSNP));
    }

    /**
     * Test of removeOffTargetVariants method, of class ExomiserSettings.
     */
    @Test
    public void testThatBuilderProducesRemoveOffTargetVariantsDefault() {
        ExomiserSettings settings = builder.build();
        assertThat(settings.removeOffTargetVariants(), is(REMOVE_OFF_TARGET_VARIANTS_DEFAULT));
    }

    @Test
    public void testThatBuilderProducesRemoveOffTargetVariantsWhenSet() {
        builder.removeOffTargetVariants(REMOVE_OFF_TARGET_VARIANTS);
        ExomiserSettings settings = builder.build();
        assertThat(settings.removeOffTargetVariants(), is(REMOVE_OFF_TARGET_VARIANTS));
    }

    /**
     * Test of getCandidateGene method, of class ExomiserSettings.
     */
    @Test
    public void testThatBuilderProducesCandidateGeneDefault() {
        ExomiserSettings settings = builder.build();
        assertThat(settings.getCandidateGene(), equalTo(CANDIDATE_GENE_NAME_DEFAULT));
    }

    @Test
    public void testThatBuilderProducesCandidateGeneWhenSet() {
        builder.candidateGene(CANDIDATE_GENE_NAME);
        ExomiserSettings settings = builder.build();
        assertThat(settings.getCandidateGene(), equalTo(CANDIDATE_GENE_NAME));
    }

    /**
     * Test of getModeOfInheritance method, of class ExomiserSettings.
     */
    @Test
    public void testThatBuilderProducesGetModeOfInheritanceDefault() {
        ExomiserSettings settings = builder.build();
        assertThat(settings.getModeOfInheritance(), equalTo(MODE_OF_INHERITANCE_DEFAULT));
    }

    @Test
    public void testThatBuilderProducesGetModeOfInheritanceWhenSet() {
        builder.modeOfInheritance(MODE_OF_INHERITANCE);
        ExomiserSettings settings = builder.build();
        assertThat(settings.getModeOfInheritance(), equalTo(MODE_OF_INHERITANCE));
    }

    /**
     * Test of getDiseaseId method, of class ExomiserSettings.
     */
    @Test
    public void testThatBuilderProducesDefaultEmptyDiseaseId() {
        ExomiserSettings settings = builder.build();
        assertThat(settings.getDiseaseId(), equalTo(DISEASE_STRING_DEFAULT));
    }

    @Test
    public void testThatBuilderProducesDiseaseIdWhenSet() {
        builder.diseaseId(DISEASE_STRING);
        ExomiserSettings settings = builder.build();
        assertThat(settings.getDiseaseId(), equalTo(DISEASE_STRING));
    }

    /**
     * Test of getHpoIds method, of class ExomiserSettings.
     */
    @Test
    public void testThatBuilderProducesDefaultEmptyHpoIds() {
        ExomiserSettings settings = builder.build();
        assertThat(settings.getHpoIds(), equalTo(HPO_LIST_DEFAULT));
    }

    @Test
    public void testThatBuilderProducesHpoIdsWhenSet() {
        builder.hpoIdList(HPO_LIST);
        ExomiserSettings settings = builder.build();
        assertThat(settings.getHpoIds(), equalTo(HPO_LIST));
    }

    /**
     * Test of getSeedGeneList method, of class ExomiserSettings.
     */
    @Test
    public void testThatBuilderProducesDefaultEmptySeedGeneList() {
        ExomiserSettings settings = builder.build();
        assertThat(settings.getSeedGeneList(), equalTo(SEED_GENE_LIST_DEFAULT));
    }

    @Test
    public void testThatBuilderProducesSeedGeneListWhenSet() {
        builder.seedGeneList(SEED_GENE_LIST);
        ExomiserSettings settings = builder.build();
        assertThat(settings.getSeedGeneList(), equalTo(SEED_GENE_LIST));
    }

    /**
     * Test of getNumberOfGenesToShow method, of class ExomiserSettings.
     */
    @Test
    public void testThatBuilderProducesDefaultNumberOfGenesToShow() {
        ExomiserSettings settings = builder.build();
        assertThat(settings.getNumberOfGenesToShow(), equalTo(NUMBER_OF_GENES_TO_SHOW_DEFAULT));
    }

    @Test
    public void testThatBuilderProducesSetNumberOfGenesToShow() {
        builder.numberOfGenesToShow(NUMBER_OF_GENES_TO_SHOW);
        ExomiserSettings settings = builder.build();
        assertThat(settings.getNumberOfGenesToShow(), equalTo(NUMBER_OF_GENES_TO_SHOW));
    }

    /**
     * Test of getOutFileName method, of class ExomiserSettings.
     */
    @Test
    public void testThatBuilderProducesDefaultOutFileName() {
        ExomiserSettings settings = builder.build();
        assertThat(settings.getOutFileName(), equalTo(OUT_FILE_NAME_DEFAULT));
    }

    @Test
    public void testThatBuilderProducesSetOutFileName() {
        builder.outFileName(OUT_FILE_NAME);
        ExomiserSettings settings = builder.build();
        assertThat(settings.getOutFileName(), equalTo(OUT_FILE_NAME));
    }

    /**
     * Test of getOutputFormat method, of class ExomiserSettings.
     */
    @Test
    public void testThatBuilderProducesDefaultOutputFormat() {
        ExomiserSettings settings = builder.build();
        assertThat(settings.getOutputFormat(), equalTo(OUTPUT_FORMAT_DEFAULT));
    }

    @Test
    public void testThatBuilderProducesSetOutputFormat() {
        builder.outputFormat(OUTPUT_FORMAT);
        ExomiserSettings settings = builder.build();
        assertThat(settings.getOutputFormat(), equalTo(OUTPUT_FORMAT));
    }

    @Test
    public void testThatBuilderCanSetAllValues() {

        builder.vcfFilePath(VCF_PATH)
                .pedFilePath(PED_PATH)
                .usePrioritiser(PriorityType.OMIM_PRIORITY)
                .maximumFrequency(MAXIMUM_FREQUENCY)
                .minimumQuality(MIMIMUM_QUALITY)
                .geneticInterval(GENETIC_INTERVAL)
                .includePathogenic(INCLUDE_PATHOGENIC)
                .removeDbSnp(REMOVE_DBSNP)
                .removeOffTargetVariants(REMOVE_OFF_TARGET_VARIANTS)
                .candidateGene(CANDIDATE_GENE_NAME)
                .modeOfInheritance(MODE_OF_INHERITANCE)
                .diseaseId(DISEASE_STRING)
                .hpoIdList(HPO_LIST)
                .seedGeneList(SEED_GENE_LIST)
                .numberOfGenesToShow(NUMBER_OF_GENES_TO_SHOW)
                .outFileName(OUT_FILE_NAME)
                .outputFormat(OUTPUT_FORMAT);

        ExomiserSettings settings = builder.build();

        assertThat(settings.getVcfPath(), equalTo(VCF_PATH));
        assertThat(settings.getPedPath(), equalTo(PED_PATH));
        assertThat(settings.getPrioritiserType(), equalTo(PriorityType.OMIM_PRIORITY));
        assertThat(settings.getMaximumFrequency(), equalTo(MAXIMUM_FREQUENCY));
        assertThat(settings.getMinimumQuality(), equalTo(MIMIMUM_QUALITY));
        assertThat(settings.getGeneticInterval(), equalTo(GENETIC_INTERVAL));
        assertThat(settings.includePathogenic(), is(INCLUDE_PATHOGENIC));
        assertThat(settings.removeDbSnp(), is(REMOVE_DBSNP));
        assertThat(settings.removeOffTargetVariants(), is(REMOVE_OFF_TARGET_VARIANTS));
        assertThat(settings.getCandidateGene(), equalTo(CANDIDATE_GENE_NAME));
        assertThat(settings.getModeOfInheritance(), equalTo(MODE_OF_INHERITANCE));
        assertThat(settings.getDiseaseId(), equalTo(DISEASE_STRING));
        assertThat(settings.getHpoIds(), equalTo(HPO_LIST));
        assertThat(settings.getSeedGeneList(), equalTo(SEED_GENE_LIST));
        assertThat(settings.getNumberOfGenesToShow(), equalTo(NUMBER_OF_GENES_TO_SHOW));
        assertThat(settings.getOutFileName(), equalTo(OUT_FILE_NAME));
        assertThat(settings.getOutputFormat(), equalTo(OUTPUT_FORMAT));

    }
    
    @Test
    public void testThatBuilderCanSetSomeValuesAndOthersRemainAsDefault() {

        builder.vcfFilePath(VCF_PATH)
                .pedFilePath(PED_PATH)
                .usePrioritiser(PriorityType.OMIM_PRIORITY)
                .maximumFrequency(MAXIMUM_FREQUENCY)
                .hpoIdList(HPO_LIST)
                .seedGeneList(SEED_GENE_LIST)
                .outFileName(OUT_FILE_NAME)
                .outputFormat(OUTPUT_FORMAT);

        ExomiserSettings settings = builder.build();

        assertThat(settings.getVcfPath(), equalTo(VCF_PATH));
        assertThat(settings.getPedPath(), equalTo(PED_PATH));
        assertThat(settings.getPrioritiserType(), equalTo(PriorityType.OMIM_PRIORITY));
        assertThat(settings.getMaximumFrequency(), equalTo(MAXIMUM_FREQUENCY));
        assertThat(settings.getMinimumQuality(), equalTo(MIMIMUM_QUALITY_DEFAULT));
        assertThat(settings.getGeneticInterval(), equalTo(GENETIC_INTERVAL_DEFAULT));
        assertThat(settings.includePathogenic(), is(INCLUDE_PATHOGENIC_DEFAULT));
        assertThat(settings.removeDbSnp(), is(REMOVE_DBSNP_DEFAULT));
        assertThat(settings.removeOffTargetVariants(), is(REMOVE_OFF_TARGET_VARIANTS_DEFAULT));
        assertThat(settings.getCandidateGene(), equalTo(CANDIDATE_GENE_NAME_DEFAULT));
        assertThat(settings.getModeOfInheritance(), equalTo(MODE_OF_INHERITANCE_DEFAULT));
        assertThat(settings.getDiseaseId(), equalTo(DISEASE_STRING_DEFAULT));
        assertThat(settings.getHpoIds(), equalTo(HPO_LIST));
        assertThat(settings.getSeedGeneList(), equalTo(SEED_GENE_LIST));
        assertThat(settings.getNumberOfGenesToShow(), equalTo(NUMBER_OF_GENES_TO_SHOW_DEFAULT));
        assertThat(settings.getOutFileName(), equalTo(OUT_FILE_NAME));
        assertThat(settings.getOutputFormat(), equalTo(OUTPUT_FORMAT));

    }

}