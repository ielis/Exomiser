package org.monarchinitiative.exomiser.core.genome.dao;

import htsjdk.tribble.readers.TabixReader;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.extension.ExtendWith;
import org.mockito.Mock;
import org.mockito.Mockito;
import org.mockito.junit.jupiter.MockitoExtension;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.model.pathogenicity.CapiceScore;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicityData;

import static org.hamcrest.CoreMatchers.equalTo;
import static org.hamcrest.MatcherAssert.assertThat;


@ExtendWith(MockitoExtension.class)
public class CapiceDaoTest {

    private CapiceDao instance;

    @Mock
    private TabixReader snvTabixReader;
    @Mock
    private TabixReader indelTabixReader;

    @BeforeEach
    void setUp() {
        TabixDataSource snvTabixDataSource = new TabixReaderAdaptor(snvTabixReader);
        TabixDataSource inDelTabixDataSource = new TabixReaderAdaptor(indelTabixReader);
        instance = new CapiceDao(snvTabixDataSource, inDelTabixDataSource);
    }

    private static VariantEvaluation variant(int chr, int pos, String ref, String alt) {
        return VariantEvaluation.builder(chr, pos, ref, alt).build();
    }

    @Test
    public void testGetPathogenicityDataSnvNoData() {
        Mockito.when(snvTabixReader.query("1:2-2")).thenReturn(MockTabixIterator.empty());
        PathogenicityData result = instance.getPathogenicityData(variant(1, 2, "A", "T"));
        assertThat(result, equalTo(PathogenicityData.empty()));
    }

    @Test
    public void testGetPathogenicityDataInsertionNoData() {
        Mockito.when(indelTabixReader.query("1:2-2")).thenReturn(MockTabixIterator.empty());
        PathogenicityData result = instance.getPathogenicityData(variant(1, 2, "C", "CA"));
        assertThat(result, equalTo(PathogenicityData.empty()));
    }

    @Test
    public void testGetPathogenicityDataInsertionSingleVariantAtPositionNoMatch() {
        Mockito.when(indelTabixReader.query("1:2-2")).thenReturn(MockTabixIterator.of(
                "1\t1\tA\tAT\t0.45",
                "1\t1\tA\tAC\t0.55"));

        PathogenicityData result = instance.getPathogenicityData(variant(1, 2, "A", "AG"));
        assertThat(result, equalTo(PathogenicityData.empty()));
    }

    @Test
    public void testGetPathogenicityDataSnvSingleVariantAtPositionOneMatch() {
        Mockito.when(snvTabixReader.query("1:2-2")).thenReturn(MockTabixIterator.of("1\t1\tA\tT\t0.234"));

        PathogenicityData result = instance.getPathogenicityData(variant(1, 2, "A", "T"));
        assertPathDataContainsCaddScore(result, 0.234f);
    }

    private void assertPathDataContainsCaddScore(PathogenicityData result, float score) {
        CapiceScore expected = CapiceScore.of(score);
        assertThat(result, equalTo(PathogenicityData.of(expected)));
    }

}