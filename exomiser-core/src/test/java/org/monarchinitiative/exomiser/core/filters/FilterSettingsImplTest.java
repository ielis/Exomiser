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
package org.monarchinitiative.exomiser.core.filters;

import de.charite.compbio.jannovar.mendel.ModeOfInheritance;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.monarchinitiative.exomiser.core.filters.FilterSettingsImpl.FilterSettingsBuilder;
import org.monarchinitiative.exomiser.core.model.GeneticInterval;

import java.util.LinkedHashSet;
import java.util.Set;

import static org.hamcrest.CoreMatchers.*;
import static org.hamcrest.MatcherAssert.assertThat;

/**
 *
 * @author Jules Jacobsen <jules.jacobsen@sanger.ac.uk>
 */
public class FilterSettingsImplTest {

    private FilterSettings instance;
    private FilterSettingsBuilder builder;

    private static final float MAX_FREQ_DEFAULT = 100.00f;
    private static final float MIN_QUAL_DEFAULT = 0;
    private static final GeneticInterval GENETIC_INTERVAL_DEFAULT = null;
    private static final boolean KEEP_NON_PATHOGENIC_VARIANTS_DEFAULT = false;
    private static final boolean REMOVE_FAILED_VARIANTS_DEFAULT = false;
    private static final boolean REMOVE_KNOWN_VARIANTS_DEFAULT = false;
    private static final boolean KEEP_OFF_TARGET_VARIANTS_DEFAULT = false;
    private static final Set<Integer> GENE_IDS_TO_KEEP_DEFAULT = new LinkedHashSet<>();
    private static final ModeOfInheritance MODE_OF_INHERITANCE_DEFAULT = ModeOfInheritance.ANY;

    @BeforeEach
    public void setUp() {
        builder = FilterSettingsImpl.builder();
    }

    @Test
    public void testDefaultValues() {
        instance = builder.build();
        assertThat(instance.getMaximumFrequency(), equalTo(MAX_FREQ_DEFAULT));
        assertThat(instance.getMinimumQuality(), equalTo(MIN_QUAL_DEFAULT));
        assertThat(instance.getGeneticInterval(), equalTo(GENETIC_INTERVAL_DEFAULT));
        assertThat(instance.keepNonPathogenicVariants(), equalTo(KEEP_NON_PATHOGENIC_VARIANTS_DEFAULT));
        assertThat(instance.removeFailedVariants(), equalTo(REMOVE_FAILED_VARIANTS_DEFAULT));
        assertThat(instance.removeKnownVariants(), equalTo(REMOVE_KNOWN_VARIANTS_DEFAULT));
        assertThat(instance.keepOffTargetVariants(), equalTo(KEEP_OFF_TARGET_VARIANTS_DEFAULT));
        assertThat(instance.getGenesToKeep(), equalTo(GENE_IDS_TO_KEEP_DEFAULT));
        assertThat(instance.getModeOfInheritance(), equalTo(MODE_OF_INHERITANCE_DEFAULT));
    }
    
    @Test
    public void testMaximumFrequency() {
        float maxFreq = 0.5f;
        instance = builder.maximumFrequency(maxFreq).build();
        assertThat(instance.getMaximumFrequency(), equalTo(maxFreq));
    }

    @Test
    public void testMinimumQuality() {
        float minQual = 10f;
        instance = builder.minimumQuality(minQual).build();
        assertThat(instance.getMinimumQuality(), equalTo(minQual));
    }

    @Test
    public void testGeneticInterval() {
        GeneticInterval interval = null;
        instance = builder.geneticInterval(interval).build();
        assertThat(instance.getGeneticInterval(), nullValue());
    }
    
    @Test
    public void testKeepNonPathogenic() {
        boolean keepNonPathogenic = false;
        instance = builder.keepNonPathogenic(keepNonPathogenic).build();
        assertThat(instance.keepNonPathogenicVariants(), equalTo(keepNonPathogenic));
    }

    @Test
    public void testRemoveFailedVariants() {
        boolean expected = true;
        instance = builder.removeFailedVariants(expected).build();
        assertThat(instance.removeFailedVariants(), equalTo(expected));
    }

    @Test
    public void testRemoveKnownVariants() {
        boolean removeKnownVariants = true;
        instance = builder.removeKnownVariants(removeKnownVariants).build();
        assertThat(instance.removeKnownVariants(), equalTo(removeKnownVariants));
    }
    
    @Test
    public void testKeepOffTargetVariants() {
        boolean keepOffTargetVariants = true;
        instance = builder.keepOffTargetVariants(keepOffTargetVariants).build();
        assertThat(instance.keepOffTargetVariants(), equalTo(keepOffTargetVariants));
    }
    
    @Test
    public void testGenesToKeep() {
        Set<String> genesToKeep = new LinkedHashSet<>();
        genesToKeep.add("1");
        genesToKeep.add("2");
        genesToKeep.add("3");
        instance = builder.genesToKeep(genesToKeep).build();
        assertThat(instance.getGenesToKeep(), equalTo(genesToKeep));        
    }
    
    @Test
    public void testModeOfInheritance() {
        ModeOfInheritance modeOfInheritance = ModeOfInheritance.AUTOSOMAL_DOMINANT;
        instance = builder.modeOfInheritance(modeOfInheritance).build();
        assertThat(instance.getModeOfInheritance(), equalTo(modeOfInheritance));
    }
    
    @Test
    public void testHashCode() {
        FilterSettings other = FilterSettingsImpl.builder().build();
        instance = builder.build();
        assertThat(instance.hashCode(), equalTo(other.hashCode()));
    }
    
    @Test
    public void testEquals() {
        FilterSettings other = FilterSettingsImpl.builder().build();
        instance = builder.build();
        assertThat(instance, equalTo(other));
    }
    
    @Test
    public void testNotEquals() {
        FilterSettings other = FilterSettingsImpl.builder().minimumQuality(Float.MAX_VALUE).build();
        instance = builder.build();
        assertThat(instance.equals(other), is(false));
    }
 
    @Test
    public void testToString() {
        instance = builder.build();
        System.out.println(instance);
    }
}
