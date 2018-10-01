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
package org.monarchinitiative.exomiser.core.prioritisers;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.hamcrest.CoreMatchers.equalTo;
import static org.hamcrest.MatcherAssert.assertThat;

/**
 *
 * @author Jules Jacobsen (jules.jacobsen@sanger.ac.uk)
 */
public class ExomeWalkerPriorityScoreTest {
    
    private ExomeWalkerPriorityResult instance;
    private final double score = 1.0;
    
    @BeforeEach
    public void setUp() {
        instance =  new ExomeWalkerPriorityResult(0, "", score);
    }

    @Test
    public void testGetPriorityType() {
        assertThat(instance.getPriorityType(), equalTo(PriorityType.EXOMEWALKER_PRIORITY));
    }

    @Test
    public void testGetScore() {
        assertThat(instance.getScore(), equalTo(score));
    }

    @Test
    public void testGetHTMLCode() {
    }

}
