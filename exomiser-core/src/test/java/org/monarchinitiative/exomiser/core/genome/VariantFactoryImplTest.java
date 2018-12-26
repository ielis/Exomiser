/*
 * The Exomiser - A tool to annotate and prioritize genomic variants
 *
 * Copyright (c) 2016-2018 Queen Mary University of London.
 * Copyright (c) 2012-2016 Charité Universitätsmedizin Berlin and Genome Research Ltd.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.monarchinitiative.exomiser.core.genome;

import com.google.common.collect.ImmutableMap;
import de.charite.compbio.jannovar.annotation.VariantEffect;
import de.charite.compbio.jannovar.data.JannovarData;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;
import org.monarchinitiative.exomiser.core.genome.jannovar.JannovarDataProtoSerialiser;
import org.monarchinitiative.exomiser.core.model.*;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.function.Consumer;
import java.util.stream.Stream;

import static java.util.stream.Collectors.toList;
import static org.hamcrest.CoreMatchers.equalTo;
import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.MatcherAssert.assertThat;

/**
 * @author Jules Jacobsen <jules.jacobsen@sanger.ac.uk>
 */
public class VariantFactoryImplTest {

    private final VariantFactory instance;

    public VariantFactoryImplTest() {
        VariantAnnotator variantAnnotator = TestFactory.buildDefaultVariantAnnotator();
        instance = new VariantFactoryImpl(variantAnnotator);
    }

    private Consumer<VariantEvaluation> printVariant() {
        return variant -> {
            GenotypesContext genotypes = variant.getVariantContext().getGenotypes();
            List<GenotypeType> genotypeTypes = genotypes.stream().map(Genotype::getType).collect(toList());
            System.out.printf("%s %s %s %s %s %s %s gene={%s %s} %s%n", variant.getChromosome(), variant.getPosition(), variant
                            .getRef(), variant.getAlt(), variant.getGenotypeString(), genotypes, genotypeTypes,
                    variant.getGeneSymbol(), variant.getGeneId(), variant.getVariantContext());
        };
    }

    @Test
    public void testStreamCreateVariantsSingleAlleles() {
        Path vcfPath = Paths.get("src/test/resources/smallTest.vcf");
        long numVariants;
        try (Stream<VariantEvaluation> variants = instance.createVariantEvaluations(vcfPath)) {
            numVariants = variants
                    .peek(printVariant())
                    .count();
        }
        assertThat(numVariants, equalTo(3L));
    }

    @Test
    public void testCreateVariantContextsMultipleAllelesDiferentSingleSampleGenotypes() {
        Path vcfPath = Paths.get("src/test/resources/multiAlleleGenotypes.vcf");
        List<VariantEvaluation> variants = instance.createVariantEvaluations(vcfPath).collect(toList());
        variants.forEach(printVariant());
        assertThat(variants.size(), equalTo(11));
    }

    @Test
    public void testCreateVariantsSingleAlleles() {
        Path vcfPath = Paths.get("src/test/resources/smallTest.vcf");
        List<VariantEvaluation> variants = instance.createVariantEvaluations(vcfPath).collect(toList());
        variants.forEach(printVariant());
        assertThat(variants.size(), equalTo(3));

    }

    @Test
    public void testCreateVariantsMultipleAllelesProduceOneVariantPerAllele() {
        Path vcfPath = Paths.get("src/test/resources/altAllele.vcf");
        List<VariantEvaluation> variants = instance.createVariantEvaluations(vcfPath).collect(toList());
        variants.forEach(printVariant());
        assertThat(variants.size(), equalTo(2));
    }

    @Test
    public void testCreateVariantsMultipleAllelesSingleSampleGenotypesShouldOnlyReturnRepresentedVariationFromGenotype() {
        Path vcfPath = Paths.get("src/test/resources/multiAlleleGenotypes.vcf");
        List<VariantEvaluation> variants = instance.createVariantEvaluations(vcfPath).collect(toList());
        variants.forEach(printVariant());
        assertThat(variants.size(), equalTo(11));
    }

    @Test
    public void testCreateVariantsNoVariantAnnotationsProduceVariantEvaluationsWithNoAnnotations() {
        Path vcfPath = Paths.get("src/test/resources/noAnnotations.vcf");
        List<VariantEvaluation> variants = instance.createVariantEvaluations(vcfPath).collect(toList());
        assertThat(variants.size(), equalTo(2));

        for (VariantEvaluation variant : variants) {
            System.out.println(variant.getChromosomeName() + " " + variant);
            assertThat(variant.hasTranscriptAnnotations(), is(false));
        }
    }

    @Test
    public void testStreamVariantEvaluationsMultipleAllelesDifferentSingleSampleGenotypes() {
        Path vcfPath = Paths.get("src/test/resources/multiAlleleGenotypes.vcf");
        List<VariantEvaluation> variants = instance.createVariantEvaluations(vcfPath).collect(toList());
        assertThat(variants.size(), equalTo(11));
    }

    @Test
    public void testKnownSingleSampleSnp() {
        Stream<VariantContext> variantContexts = TestVcfParser.forSamples("Sample")
                .parseVariantContext("10\t123256215\t.\tT\tG\t100\tPASS\tGENE=FGFR2;INHERITANCE=AD;MIM=101600\tGT\t1|0");
        List<VariantEvaluation> variants = instance.createVariantEvaluations(variantContexts)
                .peek(printVariant())
                .collect(toList());
        assertThat(variants.size(), equalTo(1));
        VariantEvaluation variantEvaluation = variants.get(0);
        System.out.println(variantEvaluation);

        assertThat(variantEvaluation.getChromosome(), equalTo(10));
        assertThat(variantEvaluation.getChromosomeName(), equalTo("10"));
        assertThat(variantEvaluation.getPosition(), equalTo(123256215));
        assertThat(variantEvaluation.getRef(), equalTo("T"));
        assertThat(variantEvaluation.getAlt(), equalTo("G"));
        assertThat(variantEvaluation.hasTranscriptAnnotations(), is(true));
        System.out.println(variantEvaluation.getTranscriptAnnotations());
        assertThat(variantEvaluation.getGeneId(), equalTo("2263"));
        assertThat(variantEvaluation.getGeneSymbol(), equalTo("FGFR2"));
        assertThat(variantEvaluation.getVariantEffect(), equalTo(VariantEffect.MISSENSE_VARIANT));
        assertThat(variantEvaluation.getSampleGenotypes(), equalTo(ImmutableMap.of("Sample", SampleGenotype.phased(AlleleCall.ALT, AlleleCall.REF))));
    }

    @Test
    public void testUnKnownSingleSampleSnp() {
        Stream<VariantContext> variantContexts = TestVcfParser.forSamples("Sample")
                .parseVariantContext("UNKNOWN\t12345\t.\tT\tC\t0\tPASS\t.\tGT:DP\t0/1:21");
        List<VariantEvaluation> variants = instance.createVariantEvaluations(variantContexts)
                .peek(printVariant())
                .collect(toList());
        assertThat(variants.size(), equalTo(1));
        VariantEvaluation variantEvaluation = variants.get(0);
        System.out.println(variantEvaluation);
        assertThat(variantEvaluation.getChromosome(), equalTo(0));
        assertThat(variantEvaluation.getChromosomeName(), equalTo("UNKNOWN"));
        assertThat(variantEvaluation.getPosition(), equalTo(12345));
        assertThat(variantEvaluation.getRef(), equalTo("T"));
        assertThat(variantEvaluation.getAlt(), equalTo("C"));
        assertThat(variantEvaluation.hasTranscriptAnnotations(), is(false));
        System.out.println(variantEvaluation.getTranscriptAnnotations());
        assertThat(variantEvaluation.getGeneId(), equalTo(""));
        assertThat(variantEvaluation.getGeneSymbol(), equalTo("."));
        assertThat(variantEvaluation.getVariantEffect(), equalTo(VariantEffect.SEQUENCE_VARIANT));
        assertThat(variantEvaluation.getSampleGenotypes(), equalTo(ImmutableMap.of("Sample", SampleGenotype.of(AlleleCall.REF, AlleleCall.ALT))));
    }

    /**
     * See https://github.com/exomiser/Exomiser/issues/207#issuecomment-310123621
     */
    @Test
    public void testSingleSampleDeletion() {
        Stream<VariantContext> variantContexts = TestVcfParser.forSamples("Sample")
                .parseVariantContext("1\t123256213\t.\tCA\tC\t100.15\tPASS\tGENE=RBM8A\tGT:DP\t0/1:33");
        List<VariantEvaluation> variants = instance.createVariantEvaluations(variantContexts)
                .peek(printVariant())
                .collect(toList());
        assertThat(variants.size(), equalTo(1));
        VariantEvaluation variantEvaluation = variants.get(0);
        System.out.println(variantEvaluation);
        assertThat(variantEvaluation.getChromosome(), equalTo(1));
        assertThat(variantEvaluation.getChromosomeName(), equalTo("1"));
        assertThat(variantEvaluation.getPosition(), equalTo(123256213));
        assertThat(variantEvaluation.getRef(), equalTo("CA"));
        assertThat(variantEvaluation.getAlt(), equalTo("C"));
        assertThat(variantEvaluation.hasTranscriptAnnotations(), is(true));
        System.out.println(variantEvaluation.getTranscriptAnnotations());
        assertThat(variantEvaluation.getGeneId(), equalTo("9939"));
        assertThat(variantEvaluation.getGeneSymbol(), equalTo("RBM8A"));
        assertThat(variantEvaluation.getVariantEffect(), equalTo(VariantEffect.INTERGENIC_VARIANT));
        assertThat(variantEvaluation.getSampleGenotypes(), equalTo(ImmutableMap.of("Sample", SampleGenotype.of(AlleleCall.REF, AlleleCall.ALT))));
    }

    @Test
    public void testSnpWithNoGenotypeReturnsNothing() {
        Stream<VariantContext> variantContexts = TestVcfParser.forSamples()
                .parseVariantContext("UNKNOWN\t12345\t.\tT\tC\t0\tPASS\t.");
        List<VariantEvaluation> variants = instance.createVariantEvaluations(variantContexts)
                .peek(printVariant())
                .collect(toList());
        assertThat(variants.isEmpty(), is(true));
    }

    @Test
    public void testSingleSampleHomVar() {
        Stream<VariantContext> variantContexts = TestVcfParser.forSamples("Sample")
                .parseVariantContext("1\t120612040\t.\tT\tTCCGCCG\t258.62\tPASS\t.\tGT\t1/1");
        List<VariantEvaluation> variants = instance.createVariantEvaluations(variantContexts)
                .peek(printVariant())
                .collect(toList());
        assertThat(variants.size(), equalTo(1));
        VariantEvaluation variant = variants.get(0);
        assertThat(variant.getSampleGenotypes(), equalTo(ImmutableMap.of("Sample", SampleGenotype.of(AlleleCall.ALT, AlleleCall.ALT))));
    }

    @Test
    public void testSingleSampleMultiPosition() {
        Stream<VariantContext> variantContexts = TestVcfParser.forSamples("Sample")
                .parseVariantContext("1\t120612040\t.\tT\tTCCGCCG,TCCTCCGCCG\t258.62\tPASS\t.\tGT\t1/2");
        List<VariantEvaluation> variants = instance.createVariantEvaluations(variantContexts)
                .peek(printVariant())
                .collect(toList());
        assertThat(variants.size(), equalTo(2));
        assertThat(variants.get(0).getSampleGenotypes(), equalTo(ImmutableMap.of("Sample", SampleGenotype.of(AlleleCall.OTHER_ALT, AlleleCall.ALT))));
        assertThat(variants.get(1).getSampleGenotypes(), equalTo(ImmutableMap.of("Sample", SampleGenotype.of(AlleleCall.OTHER_ALT, AlleleCall.ALT))));
    }

    @Test
    public void testMultiSampleMultiPositionAlleleIsSplitIntoAlternateAlleles() {
        Stream<VariantContext> variantContexts = TestVcfParser.forSamples("Sample1", "Sample2")
                .parseVariantContext("1\t120612040\t.\tT\tTCCGCCG,TCCTCCGCCG\t258.62\tPASS\t.\tGT\t0/1\t0/2");
        List<VariantEvaluation> variants = instance.createVariantEvaluations(variantContexts)
                .peek(printVariant())
                .collect(toList());
        assertThat(variants.size(), equalTo(2));
        VariantEvaluation firstAllele = variants.get(0);
        System.out.println(firstAllele);
        assertThat(firstAllele.getChromosome(), equalTo(1));
        assertThat(firstAllele.getChromosomeName(), equalTo("1"));
        assertThat(firstAllele.getPosition(), equalTo(120612040));
        assertThat(firstAllele.getRef(), equalTo("T"));
        assertThat(firstAllele.getAlt(), equalTo("TCCGCCG"));
        assertThat(firstAllele.hasTranscriptAnnotations(), is(true));
        System.out.println(firstAllele.getTranscriptAnnotations());
        assertThat(firstAllele.getGeneId(), equalTo("9939"));
        assertThat(firstAllele.getGeneSymbol(), equalTo("RBM8A"));
        assertThat(firstAllele.getVariantEffect(), equalTo(VariantEffect.INTERGENIC_VARIANT));
        assertThat(firstAllele.getSampleGenotypes(), equalTo(
                ImmutableMap.of(
                        "Sample1", SampleGenotype.of(AlleleCall.REF, AlleleCall.ALT),
                        "Sample2", SampleGenotype.of(AlleleCall.REF, AlleleCall.OTHER_ALT)
                ))
        );

        VariantEvaluation secondAllele = variants.get(1);
        System.out.println(secondAllele);
        assertThat(secondAllele.getChromosome(), equalTo(1));
        assertThat(secondAllele.getChromosomeName(), equalTo("1"));
        assertThat(secondAllele.getPosition(), equalTo(120612040));
        assertThat(secondAllele.getRef(), equalTo("T"));
        assertThat(secondAllele.getAlt(), equalTo("TCCTCCGCCG"));
        assertThat(secondAllele.hasTranscriptAnnotations(), is(true));
        System.out.println(secondAllele.getTranscriptAnnotations());
        assertThat(secondAllele.getGeneId(), equalTo("9939"));
        assertThat(secondAllele.getGeneSymbol(), equalTo("RBM8A"));
        assertThat(secondAllele.getVariantEffect(), equalTo(VariantEffect.INTERGENIC_VARIANT));
        assertThat(secondAllele.getSampleGenotypes(), equalTo(
                ImmutableMap.of(
                        "Sample1", SampleGenotype.of(AlleleCall.REF, AlleleCall.OTHER_ALT),
                        "Sample2", SampleGenotype.of(AlleleCall.REF, AlleleCall.ALT)
                ))
        );
    }

    @Test
    public void testMultiSampleMultiPositionOnlyOneAltAlleleIsPresentInSamplesProducesOneVariantEvaluation() {
        Stream<VariantContext> variantContexts = TestVcfParser.forSamples("Sample1", "Sample2")
                .parseVariantContext("1\t120612040\t.\tT\tTCCGCCG,TCCTCCGCCG\t258.62\tPASS\t.\tGT\t0/1\t0/1");
        List<VariantEvaluation> variants = instance.createVariantEvaluations(variantContexts)
                .peek(printVariant())
                .collect(toList());
        assertThat(variants.size(), equalTo(1));
        VariantEvaluation variantEvaluation = variants.get(0);
        System.out.println(variantEvaluation);
        assertThat(variantEvaluation.getChromosome(), equalTo(1));
        assertThat(variantEvaluation.getChromosomeName(), equalTo("1"));
        assertThat(variantEvaluation.getPosition(), equalTo(120612040));
        assertThat(variantEvaluation.getRef(), equalTo("T"));
        assertThat(variantEvaluation.getAlt(), equalTo("TCCGCCG"));
        assertThat(variantEvaluation.hasTranscriptAnnotations(), is(true));
        System.out.println(variantEvaluation.getTranscriptAnnotations());
        assertThat(variantEvaluation.getGeneId(), equalTo("9939"));
        assertThat(variantEvaluation.getGeneSymbol(), equalTo("RBM8A"));
        assertThat(variantEvaluation.getVariantEffect(), equalTo(VariantEffect.INTERGENIC_VARIANT));
        assertThat(variantEvaluation.getSampleGenotypes(), equalTo(
                ImmutableMap.of(
                        "Sample1", SampleGenotype.of(AlleleCall.REF, AlleleCall.ALT),
                        "Sample2", SampleGenotype.of(AlleleCall.REF, AlleleCall.ALT)
                ))
        );
    }

    @Test
    public void testStructuralVariant() {
        Stream<VariantContext> variantContexts = TestVcfParser.forSamples("Sample1")
                .parseVariantContext("10 123256215 . T <DEL> 6 PASS SVTYPE=DEL;END=123256110;SVLEN=-205;CIPOS=-56,20;CIEND=-10,62 GT:GQ 0/1:12");
        List<VariantEvaluation> variants = instance.createVariantEvaluations(variantContexts)
                .peek(printVariant())
                .collect(toList());
        assertThat(variants.size(), equalTo(1));
        VariantEvaluation variantEvaluation = variants.get(0);
        assertThat(variantEvaluation.getVariantEffect(), equalTo(VariantEffect.STRUCTURAL_VARIANT));
    }

    /**
     * Comparative performance test for loading a full genome. Ignored by default as this takes a few minutes.
     */
    @Disabled
    @Test
    public void testGenome() {

        VariantAnnotator variantAnnotator = new StubVariantAnnotator();
        VariantFactory variantFactory = new VariantFactoryImpl(variantAnnotator);

        Path vcfPath = Paths.get("C:/Users/hhx640/Documents/exomiser-cli-dev/examples/NA19722_601952_AUTOSOMAL_RECESSIVE_POMP_13_29233225_5UTR_38.vcf.gz");
        long numVariants;
        try (Stream<VariantEvaluation> variants = variantFactory.createVariantEvaluations(vcfPath)) {
            numVariants = variants
                    .count();
        }
        System.out.println("Read " + numVariants + " variants");


        VariantAnnotator jannovarVariantAnnotator = new JannovarVariantAnnotator(GenomeAssembly.HG19, loadJannovarData(), ChromosomalRegionIndex
                .empty());
        VariantFactory jannovarVariantFactory = new VariantFactoryImpl(jannovarVariantAnnotator);

        long numJannovarVariants;
        try (Stream<VariantEvaluation> variants = jannovarVariantFactory.createVariantEvaluations(vcfPath)) {
            numJannovarVariants = variants
                    .count();
        }
        System.out.println("Read " + numJannovarVariants + " variants");

    }

    private class StubVariantAnnotator implements VariantAnnotator {

        @Override
        public VariantAnnotation annotate(String chr, int pos, String ref, String alt) {
            return VariantAnnotation.builder()
                    .chromosomeName(chr)
                    .chromosome(toChromosomeNumber(chr))
                    .position(pos)
                    .ref(ref)
                    .alt(alt)
                    .build();
        }

        private int toChromosomeNumber(String chr) {
            switch (chr) {
                case "X":
                    return 23;
                case "Y":
                    return 24;
                case "M":
                case "MT":
                    return 25;
                default:
                    return Integer.parseInt(chr);
            }
        }
    }

    private JannovarData loadJannovarData() {
        Path transcriptFilePath = Paths.get("C:/Users/hhx640/Documents/exomiser-data/1806_hg19/1806_hg19_transcripts_ucsc.ser");
        return JannovarDataProtoSerialiser.load(transcriptFilePath);
    }
}
