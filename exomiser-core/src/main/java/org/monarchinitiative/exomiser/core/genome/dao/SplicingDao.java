package org.monarchinitiative.exomiser.core.genome.dao;

import de.charite.compbio.jannovar.annotation.VariantEffect;
import org.monarchinitiative.exomiser.core.model.TranscriptAnnotation;
import org.monarchinitiative.exomiser.core.model.Variant;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicityData;
import org.monarchinitiative.exomiser.core.model.pathogenicity.SplicingScore;
import org.monarchinitiative.threes.core.data.SplicingTranscriptSource;
import org.monarchinitiative.threes.core.model.GenomeCoordinates;
import org.monarchinitiative.threes.core.model.SequenceInterval;
import org.monarchinitiative.threes.core.model.SplicingTranscript;
import org.monarchinitiative.threes.core.model.SplicingVariant;
import org.monarchinitiative.threes.core.reference.fasta.GenomeSequenceAccessor;
import org.monarchinitiative.threes.core.scoring.SplicingEvaluator;
import org.monarchinitiative.threes.core.scoring.SplicingPathogenicityData;

import java.util.Optional;
import java.util.regex.Pattern;

/**
 *
 */
public class SplicingDao implements PathogenicityDao {

    /**
     * We only score variants with REF and ALT alleles that match this pattern
     */
    private static final Pattern ALLELE_REGEXP = Pattern.compile("[ACGTacgt]+");

    private final GenomeSequenceAccessor genomeSequenceAccessor;

    private final SplicingTranscriptSource splicingTranscriptSource;

    private final SplicingEvaluator splicingEvaluator;

    public SplicingDao(GenomeSequenceAccessor genomeSequenceAccessor,
                       SplicingTranscriptSource splicingTranscriptSource,
                       SplicingEvaluator splicingEvaluator) {
        this.genomeSequenceAccessor = genomeSequenceAccessor;
        this.splicingTranscriptSource = splicingTranscriptSource;
        this.splicingEvaluator = splicingEvaluator;
    }

    /**
     * Turn Exomiser's variant representation into 3S's variant object.
     *
     * @param variant Exomiser variant representation
     * @return 3S representation of the variant
     */
    private static SplicingVariant makeVariant(Variant variant) {
        return SplicingVariant.newBuilder()
                .setCoordinates(GenomeCoordinates.newBuilder()
                        .setContig(variant.getChromosomeName())
                        .setBegin(variant.getPosition() - 1) // 0-based position
                        .setEnd(variant.getPosition() + variant.getRef().length() - 1)
                        .setStrand(true)
                        .build())
                .setRef(variant.getRef())
                .setAlt(variant.getAlt())
                .build();
    }

    @Override
    public PathogenicityData getPathogenicityData(Variant variant) {
        // TODO - add check for structural variant

        // from some reason variants with '*' allele are submitted to be evaluated here
        if (!(ALLELE_REGEXP.matcher(variant.getAlt()).matches() && ALLELE_REGEXP.matcher(variant.getRef()).matches())) {
            return PathogenicityData.empty();
        }

        // we require functional annotation to be performed before calling this method, since we do not score off transcript
        // variants
        if (!variant.hasTranscriptAnnotations()) {
            return PathogenicityData.empty();
        }

        VariantEffect ve = variant.getVariantEffect();
        if (ve.isOffTranscript()) {
            return PathogenicityData.empty();
        }

        Optional<TranscriptAnnotation> txOpt = variant.getTranscriptAnnotations().stream()
                .filter(tx -> tx.getVariantEffect().equals(ve))
                .findFirst();

        if (!txOpt.isPresent()) {
            // this should not happen since variant should not have variant effect without transcript annotation in the
            // first place, however, we check anyway
            return PathogenicityData.empty();
        }

        // splicing annotations will be calculated with respect to this transcript accession
        String txAccession = txOpt.get().getAccession();

        Optional<SplicingTranscript> transcriptOptional = splicingTranscriptSource.fetchTranscripts(variant.getChromosomeName(), variant.getPosition(), variant.getPosition()).stream()
                .filter(t -> t.getAccessionId().equals(txAccession))
                .findFirst();

        if (!transcriptOptional.isPresent()) {
            // we did not find appropriate transcript record in splicing database, there's nothing we can do
            return PathogenicityData.empty();
        }

        final SplicingTranscript transcript = transcriptOptional.get();

        // get FASTA sequence of transcript
        SequenceInterval sequenceInterval = genomeSequenceAccessor.fetchSequence(variant.getChromosomeName(),
                transcript.getTxBegin() - 50, transcript.getTxEnd() + 50,
                transcript.getStrand());


        final SplicingVariant splv = makeVariant(variant);
        final SplicingPathogenicityData splicingPathogenicityData = splicingEvaluator.evaluate(splv, transcript, sequenceInterval);
        final Double maxScore = splicingPathogenicityData.getMaxScore();

        return maxScore.isNaN()
                ? PathogenicityData.empty()
                : PathogenicityData.of(SplicingScore.of(maxScore.floatValue()));
    }
}
