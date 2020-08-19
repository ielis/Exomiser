package org.monarchinitiative.exomiser.core.genome.dao;

import org.monarchinitiative.exomiser.core.model.Variant;
import org.monarchinitiative.exomiser.core.model.pathogenicity.CapiceScore;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicityData;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.cache.annotation.Cacheable;
import org.springframework.cache.annotation.Caching;

public class CapiceDao extends TabixSnvIndelDao {

    private final Logger logger = LoggerFactory.getLogger(CapiceDao.class);

    public CapiceDao(TabixDataSource capiceSnvTabixDataSource, TabixDataSource capiceInDelTabixDataSource) {
        super(capiceSnvTabixDataSource, capiceInDelTabixDataSource);
    }

    @Caching(cacheable = {
            @Cacheable(cacheNames = "hg19.capice", keyGenerator = "variantKeyGenerator", condition = "#variant.genomeAssembly == T(org.monarchinitiative.exomiser.core.genome.GenomeAssembly).HG19"),
            @Cacheable(cacheNames = "hg38.capice", keyGenerator = "variantKeyGenerator", condition = "#variant.genomeAssembly == T(org.monarchinitiative.exomiser.core.genome.GenomeAssembly).HG38"),
    })
    @Override
    public PathogenicityData getPathogenicityData(Variant variant) {
        logger.debug("Getting CAPICE data for {}", variant);
        return processResults(variant);
    }

    @Override
    protected PathogenicityData makePathogenicityData(String[] tokens) {
        /*
        CAPICE line looks like this:
        1       8390995 A       C       0.5460
        1       8390995 A       G       0.5542
        1       8390995 A       T       0.4132
         */
        float score = Float.parseFloat(tokens[4]);
        CapiceScore capiceScore = CapiceScore.of(score);
        return PathogenicityData.of(capiceScore);
    }
}
