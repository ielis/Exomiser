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

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.monarchinitiative.exomiser.core.genome.dao;

import org.monarchinitiative.exomiser.core.model.Variant;
import org.monarchinitiative.exomiser.core.model.pathogenicity.CaddScore;
import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicityData;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.cache.annotation.Cacheable;
import org.springframework.cache.annotation.Caching;

/**
 * @author Jules Jacobsen <jules.jacobsen@sanger.ac.uk>
 */
public class CaddDao extends TabixSnvIndelDao {

    private final Logger logger = LoggerFactory.getLogger(CaddDao.class);

    public CaddDao(TabixDataSource caddInDelTabixDataSource, TabixDataSource caddSnvTabixDataSource) {
        super(caddSnvTabixDataSource, caddInDelTabixDataSource);
    }

    @Caching(cacheable = {
            @Cacheable(cacheNames = "hg19.cadd", keyGenerator = "variantKeyGenerator", condition = "#variant.genomeAssembly == T(org.monarchinitiative.exomiser.core.genome.GenomeAssembly).HG19"),
            @Cacheable(cacheNames = "hg38.cadd", keyGenerator = "variantKeyGenerator", condition = "#variant.genomeAssembly == T(org.monarchinitiative.exomiser.core.genome.GenomeAssembly).HG38"),
    })
    @Override
    public PathogenicityData getPathogenicityData(Variant variant) {
        logger.debug("Getting CADD data for {}", variant);
        return processResults(variant);
    }

    @Override
    protected PathogenicityData makePathogenicityData(String[] tokens) {

        float score = Float.parseFloat(tokens[5]);
        CaddScore caddScore = CaddScore.of(score);
        return PathogenicityData.of(caddScore);
    }

}