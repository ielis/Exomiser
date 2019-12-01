package org.monarchinitiative.exomiser.core.genome.dao;

import de.charite.compbio.jannovar.annotation.VariantEffect;
import de.charite.compbio.jannovar.data.ReferenceDictionary;
import de.charite.compbio.jannovar.reference.*;
import org.monarchinitiative.exomiser.core.model.TranscriptAnnotation;
import org.monarchinitiative.exomiser.core.model.Variant;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicityData;
import org.monarchinitiative.exomiser.core.model.pathogenicity.SplicingScore;
import org.monarchinitiative.threes.core.data.SplicingTranscriptSource;
import org.monarchinitiative.threes.core.model.SplicingTranscript;
import org.monarchinitiative.threes.core.scoring.SplicingEvaluator;
import org.monarchinitiative.threes.core.scoring.SplicingPathogenicityData;
import xyz.ielis.hyperutil.reference.fasta.GenomeSequenceAccessor;
import xyz.ielis.hyperutil.reference.fasta.SequenceInterval;

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

    private final ReferenceDictionary referenceDictionary;

    public SplicingDao(GenomeSequenceAccessor genomeSequenceAccessor,
                       SplicingTranscriptSource splicingTranscriptSource,
                       SplicingEvaluator splicingEvaluator,
                       ReferenceDictionary referenceDictionary) {
        this.genomeSequenceAccessor = genomeSequenceAccessor;
        this.splicingTranscriptSource = splicingTranscriptSource;
        this.splicingEvaluator = splicingEvaluator;
        this.referenceDictionary = referenceDictionary;
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

        final Optional<SplicingTranscript> transcriptOpt = splicingTranscriptSource.fetchTranscriptByAccession(txAccession, referenceDictionary);

        if (!transcriptOpt.isPresent()) {
            // we did not find appropriate transcript record in splicing database, there's nothing we can do
            return PathogenicityData.empty();
        }

        final SplicingTranscript transcript = transcriptOpt.get();

        // get FASTA sequence of transcript
        final GenomeInterval queryInterval = transcript.getTxRegionCoordinates().withMorePadding(50, 50);
        Optional<SequenceInterval> sequenceIntervalOpt = genomeSequenceAccessor.fetchSequence(queryInterval);

        if (!sequenceIntervalOpt.isPresent()) {
            // we did not find enough FASTA sequence for given transcript
            return PathogenicityData.empty();
        }
        final SequenceInterval sequenceInterval = sequenceIntervalOpt.get();

        final GenomeVariant gv = new GenomeVariant(new GenomePosition(referenceDictionary, Strand.FWD, variant.getChromosome(), variant.getPosition(), PositionType.ONE_BASED),
                variant.getRef(), variant.getAlt());

        final SplicingPathogenicityData splicingPathogenicityData = splicingEvaluator.evaluate(gv, transcript, sequenceInterval);
        final double maxScore = splicingPathogenicityData.getMaxScore();

        return Double.isNaN(maxScore)
                ? PathogenicityData.empty()
                : PathogenicityData.of(SplicingScore.of((float) maxScore));
    }
}
