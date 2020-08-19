package org.monarchinitiative.exomiser.core.genome.dao;

import htsjdk.tribble.readers.TabixReader;
import org.monarchinitiative.exomiser.core.model.AllelePosition;
import org.monarchinitiative.exomiser.core.model.Variant;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicityData;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;

/**
 * This class operates on two TABIX files that contain genome-wide pathogenicity predictions for variants. one  one for SNVs,
 * the other for INDELs.
 */
abstract class TabixSnvIndelDao implements PathogenicityDao {
    protected final TabixDataSource snvTabixDataSource;
    protected final TabixDataSource inDelTabixDataSource;
    private final Logger logger = LoggerFactory.getLogger(CapiceDao.class);

    public TabixSnvIndelDao(TabixDataSource snvTabixDataSource, TabixDataSource inDelTabixDataSource) {
        this.snvTabixDataSource = snvTabixDataSource;
        this.inDelTabixDataSource = inDelTabixDataSource;
    }


    protected PathogenicityData processResults(Variant variant) {
        String chromosome = variant.getChromosomeName();
        String ref = variant.getRef();
        String alt = variant.getAlt();
        int start = variant.getPosition();
        if (AllelePosition.isSnv(ref, alt)) {
            return fetchPathogenicityData(snvTabixDataSource, chromosome, start, ref, alt);
        }
        return fetchPathogenicityData(inDelTabixDataSource, chromosome, start, ref, alt);
    }

    private PathogenicityData fetchPathogenicityData(TabixDataSource tabixDataSource, String chromosome, int start, String ref, String alt) {
        try {
            TabixReader.Iterator results = tabixDataSource.query(chromosome + ":" + start + "-" + start);
            String line;
            //there can be 0 - N results, we expect the TABIX files to have the following format:
            //#Chrom  Pos     Ref     Alt     ...
            //2       14962   C       CA      ...
            //2       14962   C       CAA     ...
            //2       14962   CA      C       ...
            while ((line = results.next()) != null) {
                String[] elements = line.split("\t");
                String refAllele = elements[2];
                String altAllele = elements[3];
                if (refAllele.equals(ref) && altAllele.equals(alt)) {
                    return makePathogenicityData(elements);
                }
            }
        } catch (IOException e) {
            logger.error("Unable to read from tabix file {}", tabixDataSource.getSource(), e);
        }
        return PathogenicityData.empty();
    }

    /**
     * Extract {@link PathogenicityData} from given TABIX line.
     *
     * @param tokens tokens representing the relevant TABIX line
     * @return {@link PathogenicityData} with the appropriate score
     */
    protected abstract PathogenicityData makePathogenicityData(String[] tokens);
}
