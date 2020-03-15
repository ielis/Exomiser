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

package org.monarchinitiative.exomiser.autoconfigure.genome;

import de.charite.compbio.jannovar.data.JannovarData;
import org.h2.mvstore.MVStore;
import org.monarchinitiative.exomiser.core.genome.dao.*;
import org.monarchinitiative.exomiser.core.genome.jannovar.JannovarDataSourceLoader;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.nio.file.Path;
import java.util.Objects;

/**
 * @author Jules Jacobsen <j.jacobsen@qmul.ac.uk>
 */
public class GenomeDataSourceLoader {

    private static final Logger logger = LoggerFactory.getLogger(GenomeDataSourceLoader.class);

    private final GenomeProperties genomeProperties;
    private final GenomeDataResolver genomeDataResolver;

    private final JannovarData jannovarData;
    private final MVStore mvStore;

    private final VariantWhiteList variantWhiteList;

    //TabixDataSources
    private final TabixDataSource localFrequencyTabixDataSource;
    private final TabixDataSource caddSnvTabixDataSource;
    private final TabixDataSource caddIndelTabixDataSource;
    private final TabixDataSource remmTabixDataSource;
    private final TabixDataSource testPathogenicityTabixDataSource;

    public GenomeDataSourceLoader(GenomeProperties genomeProperties, GenomeDataResolver genomeDataResolver) {
        this.genomeProperties = genomeProperties;
        this.genomeDataResolver = genomeDataResolver;

        this.jannovarData = loadJannovarData();
        this.mvStore = loadMvStore();
        this.variantWhiteList = loadVariantWhiteList();

        localFrequencyTabixDataSource = getTabixDataSourceOrDefault("LOCAL", genomeProperties.getLocalFrequencyPath());
        caddSnvTabixDataSource = getTabixDataSourceOrDefault("CADD snv", genomeProperties.getCaddSnvPath());
        caddIndelTabixDataSource = getTabixDataSourceOrDefault("CADD InDel", genomeProperties.getCaddInDelPath());
        remmTabixDataSource = getTabixDataSourceOrDefault("REMM", genomeProperties.getRemmPath());
        testPathogenicityTabixDataSource = getTabixDataSourceOrDefault("TEST", genomeProperties.getTestPathogenicityScorePath());
    }

    private JannovarData loadJannovarData() {
        Path transcriptFilePath = genomeDataResolver.getTranscriptFilePath();
        logger.debug("Loading transcript data from {}", transcriptFilePath);
        return JannovarDataSourceLoader.loadJannovarData(transcriptFilePath);
    }

    private MVStore loadMvStore() {
        Path mvStoreAbsolutePath = genomeDataResolver.getVariantsMvStorePath();
        logger.debug("Opening MVStore from {}", mvStoreAbsolutePath);
        return MvStoreDataSourceLoader.openMvStore(mvStoreAbsolutePath);
    }

    private VariantWhiteList loadVariantWhiteList() {
        Path variantWhiteListPath = genomeDataResolver.resolvePathOrNullIfEmpty(genomeProperties.getVariantWhiteListPath());
        if (variantWhiteListPath != null) {
            return VariantWhiteListLoader.loadVariantWhiteList(variantWhiteListPath);
        }
        return InMemoryVariantWhiteList.empty();
    }

    private TabixDataSource getTabixDataSourceOrDefault(String dataSourceName, String tabixPath) {
        Path path = genomeDataResolver.resolvePathOrNullIfEmpty(tabixPath);
        if (path != null) {
            logger.info("Opening {} data from source: {}", dataSourceName, path);
            return TabixDataSourceLoader.load(path);
        } else {
            logger.warn("Data for {} is not configured. THIS WILL LEAD TO ERRORS IF REQUIRED DURING ANALYSIS. Check the application.properties is pointing to a valid file.", dataSourceName);
            String message = "Data for " + dataSourceName + " is not configured. Check the application.properties is pointing to a valid file.";
            return new ErrorThrowingTabixDataSource(message);
        }
    }

    public JannovarData getJannovarData() {
        return jannovarData;
    }

    public MVStore getMvStore() {
        return mvStore;
    }

    public VariantWhiteList getVariantWhiteList() {
        return variantWhiteList;
    }

    public TabixDataSource getLocalFrequencyTabixDataSource() {
        return localFrequencyTabixDataSource;
    }

    public TabixDataSource getCaddSnvTabixDataSource() {
        return caddSnvTabixDataSource;
    }

    public TabixDataSource getCaddIndelTabixDataSource() {
        return caddIndelTabixDataSource;
    }

    public TabixDataSource getRemmTabixDataSource() {
        return remmTabixDataSource;
    }

    public TabixDataSource getTestPathogenicityTabixDataSource() {
        return testPathogenicityTabixDataSource;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        GenomeDataSourceLoader that = (GenomeDataSourceLoader) o;
        return Objects.equals(jannovarData, that.jannovarData) &&
                Objects.equals(mvStore, that.mvStore) &&
                Objects.equals(localFrequencyTabixDataSource, that.localFrequencyTabixDataSource) &&
                Objects.equals(caddSnvTabixDataSource, that.caddSnvTabixDataSource) &&
                Objects.equals(caddIndelTabixDataSource, that.caddIndelTabixDataSource) &&
                Objects.equals(remmTabixDataSource, that.remmTabixDataSource);
    }

    @Override
    public int hashCode() {
        return Objects.hash(jannovarData, mvStore, localFrequencyTabixDataSource, caddSnvTabixDataSource, caddIndelTabixDataSource, remmTabixDataSource);
    }

}
