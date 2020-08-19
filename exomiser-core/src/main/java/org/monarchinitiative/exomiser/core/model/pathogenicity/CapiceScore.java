package org.monarchinitiative.exomiser.core.model.pathogenicity;

public class CapiceScore extends BasePathogenicityScore {

    public static CapiceScore of(float score) {
        // TODO: 19. 8. 2020 is this a scaled score?
        return new CapiceScore(score);
    }

    private CapiceScore(float score) {
        super(PathogenicitySource.CAPICE, score);
    }
}
