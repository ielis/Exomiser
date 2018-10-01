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

package org.monarchinitiative.exomiser.core.analysis.util;

import com.google.common.collect.ImmutableList;
import de.charite.compbio.jannovar.mendel.ModeOfInheritance;
import org.monarchinitiative.exomiser.core.model.SampleIdentifier;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.math.BigDecimal;
import java.util.*;
import java.util.function.Predicate;

import static java.util.stream.Collectors.toList;

/**
 * Non-public helper class for finding contributing alleles.
 *
 * @author Jules Jacobsen <j.jacobsen@qmul.ac.uk>
 * @since 10.0.0
 */
class ContributingAlleleCalculator {

    private static final Logger logger = LoggerFactory.getLogger(ContributingAlleleCalculator.class);

    private final SampleIdentifier probandSampleIdentifier;
    private final CompHetAlleleCalculator compHetAlleleCalculator;

    ContributingAlleleCalculator(SampleIdentifier probandSampleIdentifier, InheritanceModeAnnotator inheritanceModeAnnotator) {
        this.probandSampleIdentifier = probandSampleIdentifier;
        this.compHetAlleleCalculator = new CompHetAlleleCalculator(inheritanceModeAnnotator);
    }

    /**
     * Calculates the total priority score for the {@code VariantEvaluation} of
     * the gene. Note that for assumed autosomal recessive variants, the mean of the worst two variants scores is
     * taken, and for other modes of inheritance, the worst (highest numerical) value is taken.
     * <p>
     * Note that we <b>cannot assume that genes have been filtered for mode of
     * inheritance before this function is called. This means that we
     * need to apply separate filtering for mode of inheritance here</b>. We also
     * need to watch out for is whether a variant is homozygous or
     * not (for autosomal recessive inheritance, these variants get counted
     * twice).
     */
    protected List<VariantEvaluation> findContributingVariantsForInheritanceMode(ModeOfInheritance modeOfInheritance, List<VariantEvaluation> variantEvaluations) {
        List<VariantEvaluation> variantsCompatibleWithMode = variantEvaluations.stream()
                //It is critical only the PASS variants are used in the scoring
                .filter(VariantEvaluation::passedFilters)
                .filter(variantEvaluation -> variantEvaluation.isCompatibleWith(modeOfInheritance))
                .collect(toList());
        //note these need to be filtered for the relevant ModeOfInheritance before being checked for the contributing variants
        if (variantsCompatibleWithMode.isEmpty()) {
            return Collections.emptyList();
        }
        switch (modeOfInheritance) {
            case AUTOSOMAL_RECESSIVE:
            case X_RECESSIVE:
                return findAutosomalRecessiveContributingVariants(modeOfInheritance, variantsCompatibleWithMode);
            default:
                return findNonAutosomalRecessiveContributingVariants(modeOfInheritance, variantsCompatibleWithMode);
        }

    }

    private List<VariantEvaluation> findAutosomalRecessiveContributingVariants(ModeOfInheritance modeOfInheritance, List<VariantEvaluation> variantEvaluations) {
        if (variantEvaluations.isEmpty()) {
            return Collections.emptyList();
        }

        Optional<CompHetPair> bestCompHetPair = compHetAlleleCalculator.findCompatibleCompHetAlleles(variantEvaluations)
                .stream()
                .map(pair -> new CompHetPair(pair.get(0), pair.get(1)))
                .max(Comparator.comparing(CompHetPair::getScore));

        Optional<VariantEvaluation> bestHomozygousAlt = variantEvaluations.stream()
                .filter(variantIsHomozygousAlt(probandSampleIdentifier))
                .max(Comparator.comparing(VariantEvaluation::getVariantScore));

        // Realised original logic allows a comphet to be calculated between a top scoring het and second place hom which is wrong
        // Jannovar seems to currently be allowing hom_ref variants through so skip these as well
        double bestCompHetScore = bestCompHetPair
                .map(CompHetPair::getScore)
                .orElse(0.0);

        double bestHomAltScore = bestHomozygousAlt
                .map(VariantEvaluation::getVariantScore)
                .orElse(0f);

        double bestScore = Double.max(bestHomAltScore, bestCompHetScore);

        if (BigDecimal.valueOf(bestScore).equals(BigDecimal.valueOf(bestCompHetScore)) && bestCompHetPair.isPresent()) {
            CompHetPair compHetPair = bestCompHetPair.get();
            compHetPair.setContributesToGeneScoreUnderMode(modeOfInheritance);
            logger.debug("Top scoring comp het: {}", compHetPair);

            return bestCompHetPair.get().getAlleles();
        } else if (bestHomozygousAlt.isPresent()) {
            VariantEvaluation topHomAlt = bestHomozygousAlt.get();
            topHomAlt.setContributesToGeneScoreUnderMode(modeOfInheritance);
            logger.debug("Top scoring hom alt het: {}", topHomAlt);

            return ImmutableList.of(topHomAlt);
        }
        return Collections.emptyList();
    }

    private List<VariantEvaluation> findNonAutosomalRecessiveContributingVariants(ModeOfInheritance modeOfInheritance, List<VariantEvaluation> variantEvaluations) {
        Optional<VariantEvaluation> bestVariant = variantEvaluations.stream()
                .max(Comparator.comparing(VariantEvaluation::getVariantScore));

        bestVariant.ifPresent(variantEvaluation -> variantEvaluation.setContributesToGeneScoreUnderMode(modeOfInheritance));

        return bestVariant.map(Collections::singletonList).orElseGet(Collections::emptyList);
    }

    private Predicate<VariantEvaluation> variantIsHomozygousAlt(SampleIdentifier probandSampleIdentifier) {
        return ve -> ve.getSampleGenotype(probandSampleIdentifier.getId()).isHomozygousAlt();
    }

    /**
     * Data class for holding pairs of alleles which are compatible with AR compound heterozygous inheritance.
     */
    private static final class CompHetPair {

        private final double score;
        private final VariantEvaluation allele1;
        private final VariantEvaluation allele2;

        CompHetPair(VariantEvaluation allele1, VariantEvaluation allele2) {
            this.allele1 = allele1;
            this.allele2 = allele2;
            this.score = calculateScore(allele1, allele2);
        }

        double calculateScore(VariantEvaluation allele1, VariantEvaluation allele2) {
            double allele1Score = allele1 == null ? 0 : allele1.getVariantScore();
            double allele2Score = allele2 == null ? 0 : allele2.getVariantScore();
            return (allele1Score + allele2Score) / 2.0;
        }

        double getScore() {
            return score;
        }

        List<VariantEvaluation> getAlleles() {
            List<VariantEvaluation> alleles = new ArrayList<>();
            if (null != allele1) {
                alleles.add(allele1);
            }
            if (null != allele2) {
                alleles.add(allele2);
            }
            return alleles;
        }

        void setContributesToGeneScoreUnderMode(ModeOfInheritance modeOfInheritance) {
            if (allele1 != null) {
                allele1.setContributesToGeneScoreUnderMode(modeOfInheritance);
            }
            if (allele2 != null) {
                allele2.setContributesToGeneScoreUnderMode(modeOfInheritance);
            }
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            CompHetPair that = (CompHetPair) o;
            return Double.compare(that.score, score) == 0 &&
                    Objects.equals(allele1, that.allele1) &&
                    Objects.equals(allele2, that.allele2);
        }

        @Override
        public int hashCode() {
            return Objects.hash(allele1, allele2, score);
        }

        @Override
        public String toString() {
            return "CompHetPair{" +
                    "score=" + score +
                    ", allele1=" + allele1 +
                    ", allele2=" + allele2 +
                    '}';
        }
    }
}
