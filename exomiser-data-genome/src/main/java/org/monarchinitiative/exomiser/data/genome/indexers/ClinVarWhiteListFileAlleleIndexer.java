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

package org.monarchinitiative.exomiser.data.genome.indexers;

import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
import org.monarchinitiative.exomiser.data.genome.model.Allele;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.StringJoiner;
import java.util.concurrent.atomic.AtomicLong;

/**
 * Specialised AlleleIndexer for producing the ClinVar whitelist
 *
 * @author Jules Jacobsen <j.jacobsen@qmul.ac.uk>
 */
public class ClinVarWhiteListFileAlleleIndexer extends AbstractAlleleIndexer {

    private static final Logger logger = LoggerFactory.getLogger(ClinVarWhiteListFileAlleleIndexer.class);

    private final BufferedWriter bufferedWriter;
    private final AtomicLong count = new AtomicLong(0);

    public ClinVarWhiteListFileAlleleIndexer(BufferedWriter bufferedWriter) {
        this.bufferedWriter = bufferedWriter;
    }

    @Override
    public void writeAllele(Allele allele) {
        ClinVarData clinVarData = allele.getClinVarData();
        if (isPathOrLikelyPath(clinVarData) && hasAssertionCriteria(clinVarData)) {
            StringJoiner stringJoiner = new StringJoiner("\t");
            stringJoiner.add(Integer.toString(allele.getChr()));
            stringJoiner.add(Integer.toString(allele.getPos()));
            stringJoiner.add(allele.getRef());
            stringJoiner.add(allele.getAlt());
            stringJoiner.merge(createClinVarInfo(clinVarData));

            try {
                bufferedWriter.write(stringJoiner.toString());
                bufferedWriter.newLine();
                // when writing out to the BlockCompressedOutputStream don't use flush()
                // this is necessary if only writing to an uncompressed test file.
            } catch (IOException ex) {
                logger.error("Unable to write to allele index file", ex);
                throw new RuntimeException(ex);
            }
            count.incrementAndGet();
        }
    }

    private boolean hasAssertionCriteria(ClinVarData clinVarData) {
        // maps to the CLNREVSTAT subfield in the VCF INFO. Many alleles with 'no_assertion_criteria_provided'
        // or 'no_assertion_provided' have incredibly high MAF, some even as high as 98% in some populations.
        return !clinVarData.getReviewStatus().startsWith("no_assertion");
    }

    private boolean isPathOrLikelyPath(ClinVarData clinVarData) {
        switch(clinVarData.getPrimaryInterpretation()) {
            case PATHOGENIC:
            case PATHOGENIC_OR_LIKELY_PATHOGENIC:
            case LIKELY_PATHOGENIC:
                return true;
            default:
                return false;
        }
    }

    private StringJoiner createClinVarInfo(ClinVarData clinVarData) {
        StringJoiner stringJoiner = new StringJoiner(";");
        stringJoiner.add("ALLELEID=" + clinVarData.getAlleleId());
        stringJoiner.add("CLNSIG=" + clinVarData.getPrimaryInterpretation());
        stringJoiner.add("CLNREVSTAT=" + clinVarData.getReviewStatus());
        return stringJoiner;
    }

    @Override
    public long count() {
        return count.get();
    }

    @Override
    public void close() {
        try {
            bufferedWriter.close();
        } catch (IOException e) {
            logger.error("Unable to close allele index file", e);
            throw new RuntimeException(e);
        }
    }
}
