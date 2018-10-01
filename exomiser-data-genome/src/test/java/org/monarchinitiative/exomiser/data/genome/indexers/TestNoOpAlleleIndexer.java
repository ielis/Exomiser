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

package org.monarchinitiative.exomiser.data.genome.indexers;

import org.monarchinitiative.exomiser.data.genome.model.Allele;

/**
 * Test no-op {@link AlleleIndexer} which does nothing other than count the number of alleles it was asked to write.
 * Useful for integration testing and starting to write new {@link AlleleIndexer} classes.
 *
 * @author Jules Jacobsen <j.jacobsen@qmul.ac.uk>
 */
public class TestNoOpAlleleIndexer extends AbstractAlleleIndexer {

    private long count;

    @Override
    protected void writeAllele(Allele allele) {
        count++;
    }

    @Override
    public long count() {
        return count;
    }

    @Override
    public void close() {
        //nothing to close
    }
}
