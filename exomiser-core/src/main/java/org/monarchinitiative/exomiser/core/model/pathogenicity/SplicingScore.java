package org.monarchinitiative.exomiser.core.model.pathogenicity;

/**
 *
 */
public class SplicingScore extends BasePathogenicityScore {

    /**
     * @param scaledScore the raw score, scaled to fit the 0 - 1 scaling
     */
    private SplicingScore(float scaledScore) {
        super(PathogenicitySource.SPLICING, scaledScore);
    }

    public static SplicingScore of(float score) {
        return new SplicingScore(score);
    }

}
