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

package org.monarchinitiative.exomiser.core.analysis;

import de.charite.compbio.jannovar.mendel.ModeOfInheritance;
import org.junit.jupiter.api.Test;
import org.monarchinitiative.exomiser.core.filters.InheritanceFilter;
import org.monarchinitiative.exomiser.core.filters.KnownVariantFilter;
import org.monarchinitiative.exomiser.core.filters.PriorityScoreFilter;
import org.monarchinitiative.exomiser.core.prioritisers.OmimPriority;
import org.monarchinitiative.exomiser.core.prioritisers.PhivePriority;
import org.monarchinitiative.exomiser.core.prioritisers.PriorityType;
import org.monarchinitiative.exomiser.core.prioritisers.service.TestPriorityServiceFactory;

import java.util.EnumSet;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.core.Is.is;

/**
 * @author Jules Jacobsen <jules.jacobsen@sanger.ac.uk>
 */
public class AnalysisStepTest {

    private static final AnalysisStep PHIVE_PRIORITY = new PhivePriority(TestPriorityServiceFactory.stubPriorityService());
    private static final AnalysisStep PRIORITY_SCORE_FILTER = new PriorityScoreFilter(PriorityType.PHIVE_PRIORITY, 0.6f);
    private static final AnalysisStep KNOWN_VARIANT_FILTER = new KnownVariantFilter();
    private static final AnalysisStep OMIM_PRIORITY = new OmimPriority(TestPriorityServiceFactory.stubPriorityService());
    private static final AnalysisStep INHERITANCE_FILTER = new InheritanceFilter(EnumSet.of(ModeOfInheritance.ANY));
    
    @Test
    public void testIsInheritanceModeDependentOMIMPriority() {
        assertThat(OMIM_PRIORITY.isInheritanceModeDependent(), is(true));
    }

    @Test
    public void testIsInheritanceModeDependentInheritanceModeFilter() {
        assertThat(INHERITANCE_FILTER.isInheritanceModeDependent(), is(true));
    }

    @Test
    public void testIsInheritanceModeDependentNotInheritanceModeDependant() {
        assertThat(KNOWN_VARIANT_FILTER.isInheritanceModeDependent(), is(false));
    }

    @Test
    public void testIsOnlyGeneDependentInheritanceModeDependantGeneFilter() {
        assertThat(INHERITANCE_FILTER.isOnlyGeneDependent(), is(false));
    }

    @Test
    public void testIsOnlyGeneDependentInheritanceModeDependantPrioritiser() {
        assertThat(OMIM_PRIORITY.isOnlyGeneDependent(), is(false));
    }

    @Test
    public void testIsOnlyGeneDependentVariantFilter() {
        assertThat(KNOWN_VARIANT_FILTER.isOnlyGeneDependent(), is(false));
    }

    @Test
    public void testIsOnlyGeneDependentOtherPrioritiser() {
        assertThat(PHIVE_PRIORITY.isOnlyGeneDependent(), is(true));
    }

    @Test
    public void testIsOnlyGeneDependentPriorityScoreFilter() {
        assertThat(PRIORITY_SCORE_FILTER.isOnlyGeneDependent(), is(true));
    }

    @Test
    public void testIsVariantFilterVariantFilter() {
        assertThat(KNOWN_VARIANT_FILTER.isVariantFilter(), is(true));
    }

    @Test
    public void testIsVariantFilterGeneFilter() {
        assertThat(INHERITANCE_FILTER.isVariantFilter(), is(false));
    }

    @Test
    public void testIsVariantFilterPrioritiser() {
        assertThat(OMIM_PRIORITY.isVariantFilter(), is(false));
    }

}