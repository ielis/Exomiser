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

package org.monarchinitiative.exomiser.core.genome.dao;

import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.extension.ExtendWith;
import org.monarchinitiative.exomiser.core.model.StructuralType;
import org.monarchinitiative.exomiser.core.model.Variant;
import org.monarchinitiative.exomiser.core.model.VariantAnnotation;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicityData;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicityScore;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicitySource;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit.jupiter.SpringExtension;

import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.MatcherAssert.assertThat;

/**
 * @author Jules Jacobsen <j.jacobsen@qmul.ac.uk>
 */
@Disabled
@ExtendWith(SpringExtension.class)
@ContextConfiguration(classes = {TestSvDataSourceConfig.class, SvPathogenicityDao.class})
class SvPathogenicityDaoTest {

    @Autowired
    private SvPathogenicityDao instance;

    @Test
    void getDupExactMatch() {
        Variant variant = VariantAnnotation.builder()
                .chromosome(17)
                .start(526)
                .end(81_041_938)
                .structuralType(StructuralType.DUP)
                .build();

        PathogenicityData result = instance.getPathogenicityData(variant);

        System.out.println(result);
    }

    @Test
    void getDupInexactMatch() {
        Variant variant = VariantAnnotation.builder()
                .chromosome(17)
                .start(500)
                .end(82_041_938)
                .structuralType(StructuralType.DUP)
                .build();

        PathogenicityData result = instance.getPathogenicityData(variant);

        System.out.println(result);
    }

    @Test
    void getInsExactMatch() {
        Variant variant = VariantAnnotation.builder()
                .chromosome(10)
                .start(105_817_214)
                .end(105_817_214)
                .structuralType(StructuralType.DUP)
                .build();

        PathogenicityData result = instance.getPathogenicityData(variant);

        System.out.println(result);
    }

    @Test
    void getPathogenicityData() {
        Variant variant = VariantAnnotation.builder()
                .chromosome(9)
                .start(133_388_590)
                .end(141_018_980)
                .structuralType(StructuralType.DUP)
                .build();

        PathogenicityData result = instance.getPathogenicityData(variant);

        assertThat(result, is(PathogenicityData.of(ClinVarData.builder()
                        .alleleId("RCV000133778.3")
                        .primaryInterpretation(ClinVarData.ClinSig.PATHOGENIC)
                        .reviewStatus("")
                        .build(),
                PathogenicityScore.of(PathogenicitySource.DBVAR, 1.000F),
                PathogenicityScore.of(PathogenicitySource.ISCA, 1.000F))));
    }
}