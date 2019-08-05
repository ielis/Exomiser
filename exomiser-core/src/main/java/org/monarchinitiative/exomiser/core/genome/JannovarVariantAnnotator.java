/*
 * The Exomiser - A tool to annotate and prioritize genomic variants
 *
 * Copyright (c) 2016-2019 Queen Mary University of London.
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

package org.monarchinitiative.exomiser.core.genome;

import com.google.common.collect.ImmutableList;
import de.charite.compbio.jannovar.annotation.Annotation;
import de.charite.compbio.jannovar.annotation.PutativeImpact;
import de.charite.compbio.jannovar.annotation.VariantAnnotations;
import de.charite.compbio.jannovar.annotation.VariantEffect;
import de.charite.compbio.jannovar.data.JannovarData;
import de.charite.compbio.jannovar.hgvs.AminoAcidCode;
import de.charite.compbio.jannovar.reference.GenomeVariant;
import de.charite.compbio.jannovar.reference.TranscriptModel;
import org.monarchinitiative.exomiser.core.model.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.toList;

/**
 * Handles creation of {@link VariantAnnotation} using Jannovar.
 *
 * @author Jules Jacobsen <j.jacobsen@qmul.ac.uk>
 */
public class JannovarVariantAnnotator implements VariantAnnotator {

    private static final Logger logger = LoggerFactory.getLogger(JannovarVariantAnnotator.class);

    private final GenomeAssembly genomeAssembly;
    private final JannovarAnnotationService jannovarAnnotationService;
    private final ChromosomalRegionIndex<RegulatoryFeature> regulatoryRegionIndex;

    public JannovarVariantAnnotator(GenomeAssembly genomeAssembly, JannovarData jannovarData, ChromosomalRegionIndex<RegulatoryFeature> regulatoryRegionIndex) {
        this.genomeAssembly = genomeAssembly;
        this.jannovarAnnotationService = new JannovarAnnotationService(jannovarData);
        this.regulatoryRegionIndex = regulatoryRegionIndex;
    }

    /**
     * Given a single allele from a multi-positional site, incoming variants might not be fully trimmed.
     * In cases where there is repetition, depending on the program used, the final variant allele will be different.
     * VCF:      X-118887583-TCAAAA-TCAAAACAAAA
     * Exomiser: X-118887583-T     -TCAAAA
     * CellBase: X-118887584--     - CAAAA
     * Jannovar: X-118887588-      -      CAAAA
     * Nirvana:  X-118887589-      -      CAAAA
     * <p>
     * Trimming first with Exomiser, then annotating with Jannovar, constrains the Jannovar annotation to the same
     * position as Exomiser.
     * VCF:      X-118887583-TCAAAA-TCAAAACAAAA
     * Exomiser: X-118887583-T     -TCAAAA
     * CellBase: X-118887584--     - CAAAA
     * Jannovar: X-118887583-      - CAAAA      (Jannovar is zero-based)
     * Nirvana:  X-118887584-      - CAAAA
     * <p>
     * Cellbase:
     * https://github.com/opencb/biodata/blob/develop/biodata-tools/src/main/java/org/opencb/biodata/tools/variant/VariantNormalizer.java
     * http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/v4/hsapiens/genomic/variant/X:118887583:TCAAAA:TCAAAACAAAA/annotation?assembly=grch37&limit=-1&skip=-1&count=false&Output format=json&normalize=true
     * <p>
     * Nirvana style trimming:
     * https://github.com/Illumina/Nirvana/blob/master/VariantAnnotation/Algorithms/BiDirectionalTrimmer.cs
     * <p>
     * Jannovar:
     * https://github.com/charite/jannovar/blob/master/jannovar-core/src/main/java/de/charite/compbio/jannovar/reference/VariantDataCorrector.java
     *
     * @param chr chromosome identifier
     * @param pos 1-based start position of the first base of the ref string
     * @param ref reference base(s)
     * @param alt alternate bases
     * @return {@link VariantAnnotation} objects trimmed according to {@link AllelePosition#trim(int, String, String)} and annotated using Jannovar.
     * @since 13.0.0
     */
    @Override
    public List<VariantAnnotation> annotate(String chr, int pos, String ref, String alt) {
        //so given the above, trim the allele first, then annotate it otherwise untrimmed alleles from multi-allelic sites will give different results
        AllelePosition trimmedAllele = AllelePosition.trim(pos, ref, alt);
        VariantAnnotations variantAnnotations = jannovarAnnotationService
                .annotateVariant(chr, trimmedAllele.getPos(), trimmedAllele.getRef(), trimmedAllele.getAlt());

        // Group annotations by geneSymbol then create new Jannovar.VariantAnnotations from these then return List<VariantAnnotation>
        // see issue https://github.com/exomiser/Exomiser/issues/294. However it creates approximately 2x as many variants
        // which doubles the runtime, and most of the new variants are then filtered out. So here we're trying to limit the amount of new
        // VariantAnnotations returned by only splitting those with a MODERATE or greater putative impact.
        if (effectsMoreThanOneGeneWithMinimumImpact(variantAnnotations, PutativeImpact.MODERATE)) {
            return splitAnnotationsByGene(variantAnnotations)
                    .map(variantGeneAnnotations -> buildVariantAlleleAnnotation(genomeAssembly, chr, trimmedAllele, variantGeneAnnotations))
                    .collect(toList());
        }
        return ImmutableList.of(buildVariantAlleleAnnotation(genomeAssembly, chr, trimmedAllele, variantAnnotations));
    }

    private boolean effectsMoreThanOneGeneWithMinimumImpact(VariantAnnotations variantAnnotations, PutativeImpact minimumImpact) {
        ImmutableList<Annotation> annotations = variantAnnotations.getAnnotations();
        // increasing cost of computation with each stage of the evaluation - don't change order.
        return annotations.size() > 1 &&
                variantEffectImpactIsAtLeast(variantAnnotations.getHighestImpactEffect(), minimumImpact) &&
                annotationsContainMoreThanOneGene(annotations) &&
                annotationsEffectMoreThanOneGeneWithMinimumImpact(annotations, minimumImpact);
    }

    private boolean variantEffectImpactIsAtLeast(VariantEffect variantEffect, PutativeImpact minimumImpact) {
        // HIGH and MODERATE are ordinal 0, 1 which we want to look at if there are any
        return variantEffect.getImpact().ordinal() <= minimumImpact.ordinal();
    }

    private boolean annotationsContainMoreThanOneGene(ImmutableList<Annotation> annotations) {
        return annotations.stream().map(Annotation::getGeneSymbol).distinct().count() > 1L;
    }

    private boolean annotationsEffectMoreThanOneGeneWithMinimumImpact(ImmutableList<Annotation> annotations, PutativeImpact minimumImpact) {
        Map<String, List<Annotation>> annotationsByGene = annotations.stream()
                .filter(annotation -> annotation.getMostPathogenicVarType() != null)
                .filter(annotation -> variantEffectImpactIsAtLeast(annotation.getMostPathogenicVarType(), minimumImpact))
                .collect(groupingBy(Annotation::getGeneSymbol));
        return annotationsByGene.size() > 1;
    }

    private Stream<VariantAnnotations> splitAnnotationsByGene(VariantAnnotations variantAnnotations) {
        ImmutableList<Annotation> annotations = variantAnnotations.getAnnotations();
        GenomeVariant genomeVariant = variantAnnotations.getGenomeVariant();
        logger.debug("Multiple annotations for {} {} {} {}", genomeVariant.getChrName(), genomeVariant.getPos(), genomeVariant.getRef(), genomeVariant.getAlt());

        return annotations.stream()
                .collect(groupingBy(Annotation::getGeneSymbol))
                .values().stream()
                //.peek(annotationList -> annotationList.forEach(annotation -> logger.info("{}", toAnnotationString(annotation))))
                .map(annos -> new VariantAnnotations(genomeVariant, annos));
    }

    private String toAnnotationString(Annotation annotation) {
        return annotation.getGeneSymbol() + ", " + annotation.getMostPathogenicVarType() + ", " + annotation.getMostPathogenicVarType()
                .getImpact() + ", " + annotation.getTranscript();
    }

    private VariantAnnotation buildVariantAlleleAnnotation(GenomeAssembly genomeAssembly, String contig, AllelePosition allelePosition, VariantAnnotations variantAnnotations) {
        int chr = variantAnnotations.getChr();
        GenomeVariant genomeVariant = variantAnnotations.getGenomeVariant();
        String chromosomeName = genomeVariant.getChrName() == null ? contig : genomeVariant.getChrName();
        //Attention! highestImpactAnnotation can be null
        Annotation highestImpactAnnotation = variantAnnotations.getHighestImpactAnnotation();
        String geneSymbol = buildGeneSymbol(highestImpactAnnotation);
        String geneId = buildGeneId(highestImpactAnnotation);

        //Jannovar presently ignores all structural variants, so flag it here. Not that we do anything with them at present.
        VariantEffect highestImpactEffect = allelePosition.isSymbolic() ? VariantEffect.STRUCTURAL_VARIANT : variantAnnotations.getHighestImpactEffect();
        List<TranscriptAnnotation> annotations = buildTranscriptAnnotations(variantAnnotations.getAnnotations());

        int pos = allelePosition.getPos();
        String ref = allelePosition.getRef();
        String alt = allelePosition.getAlt();

        VariantEffect variantEffect = checkRegulatoryRegionVariantEffect(highestImpactEffect, chr, pos);
        return VariantAnnotation.builder()
                .genomeAssembly(genomeAssembly)
                .chromosome(chr)
                .chromosomeName(chromosomeName)
                .position(pos)
                .ref(ref)
                .alt(alt)
                .geneId(geneId)
                .geneSymbol(geneSymbol)
                .variantEffect(variantEffect)
                .annotations(annotations)
                .build();
    }

    private List<TranscriptAnnotation> buildTranscriptAnnotations(List<Annotation> annotations) {
        List<TranscriptAnnotation> transcriptAnnotations = new ArrayList<>(annotations.size());
        for (Annotation annotation : annotations) {
            transcriptAnnotations.add(toTranscriptAnnotation(annotation));
        }
        return transcriptAnnotations;
    }

    private TranscriptAnnotation toTranscriptAnnotation(Annotation annotation) {
        return TranscriptAnnotation.builder()
                .variantEffect(annotation.getMostPathogenicVarType())
                .accession(getTranscriptAccession(annotation))
                .geneSymbol(buildGeneSymbol(annotation))
                .hgvsGenomic((annotation.getGenomicNTChange() == null) ? "" : annotation.getGenomicNTChangeStr())
                .hgvsCdna(annotation.getCDSNTChangeStr())
                .hgvsProtein(annotation.getProteinChangeStr(AminoAcidCode.THREE_LETTER))
                .distanceFromNearestGene(getDistFromNearestGene(annotation))
                .build();
    }

    private String getTranscriptAccession(Annotation annotation) {
        TranscriptModel transcriptModel = annotation.getTranscript();
        if (transcriptModel == null) {
            return "";
        }
        return transcriptModel.getAccession();
    }

    private int getDistFromNearestGene(Annotation annotation) {

        TranscriptModel tm = annotation.getTranscript();
        if (tm == null) {
            return Integer.MIN_VALUE;
        }
        GenomeVariant change = annotation.getGenomeVariant();
        Set<VariantEffect> effects = annotation.getEffects();
        if (effects.contains(VariantEffect.INTERGENIC_VARIANT) || effects.contains(VariantEffect.UPSTREAM_GENE_VARIANT) || effects
                .contains(VariantEffect.DOWNSTREAM_GENE_VARIANT)) {
            if (change.getGenomeInterval().isLeftOf(tm.getTXRegion().getGenomeBeginPos()))
                return tm.getTXRegion().getGenomeBeginPos().differenceTo(change.getGenomeInterval().getGenomeEndPos());
            else
                return change.getGenomeInterval().getGenomeBeginPos().differenceTo(tm.getTXRegion().getGenomeEndPos());
        }

        return Integer.MIN_VALUE;
    }

    private String buildGeneId(Annotation annotation) {
        if (annotation == null) {
            return "";
        }

        final TranscriptModel transcriptModel = annotation.getTranscript();
        if (transcriptModel == null) {
            return "";
        }
        //this will now return the id from the user-specified data source. Previously would only return the Entrez id.
        String geneId = transcriptModel.getGeneID();
        return geneId == null ? "" : geneId;
    }

    private String buildGeneSymbol(Annotation annotation) {
        if (annotation == null || annotation.getGeneSymbol() == null) {
            return ".";
        } else {
            return annotation.getGeneSymbol();
        }
    }

    //Adds the missing REGULATORY_REGION_VARIANT effect to variants - this isn't in the Jannovar data set.
    private VariantEffect checkRegulatoryRegionVariantEffect(VariantEffect variantEffect, int chr, int pos) {
        //n.b this check here is important as ENSEMBLE can have regulatory regions overlapping with missense variants.
        if (isIntergenicOrUpstreamOfGene(variantEffect) && regulatoryRegionIndex.hasRegionContainingPosition(chr, pos)) {
            //the effect is the same for all regulatory regions, so for the sake of speed, just assign it here rather than look it up from the list
            return VariantEffect.REGULATORY_REGION_VARIANT;
        }
        return variantEffect;
    }

    private boolean isIntergenicOrUpstreamOfGene(VariantEffect variantEffect) {
        return variantEffect == VariantEffect.INTERGENIC_VARIANT || variantEffect == VariantEffect.UPSTREAM_GENE_VARIANT;
    }

}
