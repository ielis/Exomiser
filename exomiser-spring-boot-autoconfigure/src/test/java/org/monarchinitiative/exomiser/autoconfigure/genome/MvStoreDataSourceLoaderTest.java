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

package org.monarchinitiative.exomiser.autoconfigure.genome;

import org.h2.mvstore.MVStore;
import org.junit.jupiter.api.Test;

import java.nio.file.Path;
import java.nio.file.Paths;

import static org.hamcrest.CoreMatchers.instanceOf;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * @author Jules Jacobsen <j.jacobsen@qmul.ac.uk>
 */
public class MvStoreDataSourceLoaderTest {

    @Test
    public void loadsData() {
        Path mvStorePath = Paths.get("src/test/resources/data/1710_hg19/1710_hg19_variants.mv.db");
        MVStore mvStore = MvStoreDataSourceLoader.openMvStore(mvStorePath);
        assertThat(mvStore, instanceOf(MVStore.class));
    }

    @Test
    public void cannotLoadData() {
        Path mvStorePath = Paths.get("wibble");
        assertThrows(IllegalStateException.class, () -> MvStoreDataSourceLoader.openMvStore(mvStorePath));
    }
}