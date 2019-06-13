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

import de.charite.compbio.jannovar.data.JannovarData;
import org.h2.mvstore.MVStore;
import org.monarchinitiative.exomiser.autoconfigure.DataDirectoryAutoConfiguration;
import org.monarchinitiative.exomiser.core.genome.*;
import org.monarchinitiative.exomiser.core.genome.dao.*;
import org.monarchinitiative.threes.core.data.ContigLengthDao;
import org.monarchinitiative.threes.core.data.DbSplicingTranscriptSource;
import org.monarchinitiative.threes.core.data.SplicingTranscriptSource;
import org.monarchinitiative.threes.core.pwm.DbSplicingPositionalWeightMatrixParser;
import org.monarchinitiative.threes.core.pwm.SplicingInformationContentAnnotator;
import org.monarchinitiative.threes.core.pwm.SplicingPositionalWeightMatrixParser;
import org.monarchinitiative.threes.core.reference.GenomeCoordinatesFlipper;
import org.monarchinitiative.threes.core.reference.transcript.NaiveSplicingTranscriptLocator;
import org.monarchinitiative.threes.core.reference.transcript.SplicingTranscriptLocator;
import org.monarchinitiative.threes.core.scoring.SimpleSplicingEvaluator;
import org.monarchinitiative.threes.core.scoring.SplicingEvaluator;
import org.monarchinitiative.threes.core.scoring.scorers.RawScoringFactory;
import org.monarchinitiative.threes.core.scoring.scorers.ScorerFactory;
import org.springframework.boot.autoconfigure.condition.ConditionalOnProperty;
import org.springframework.boot.context.properties.EnableConfigurationProperties;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.context.annotation.Import;

import java.nio.file.Path;

/**
 * @author Jules Jacobsen <j.jacobsen@qmul.ac.uk>
 */
@Configuration
@Import({DataDirectoryAutoConfiguration.class})
@ConditionalOnProperty({"exomiser.hg38.data-version"})
@EnableConfigurationProperties(Hg38GenomeProperties.class)
public class Hg38GenomeAnalysisServiceAutoConfiguration extends GenomeAnalysisServiceConfigurer {

    public Hg38GenomeAnalysisServiceAutoConfiguration(Hg38GenomeProperties hg38GenomeProperties, Path exomiserDataDirectory) {
        super(hg38GenomeProperties, exomiserDataDirectory);
    }

    @Bean("hg38jannovarData")
    public JannovarData jannovarData() {
        return jannovarData;
    }

    @Bean("hg38mvStore")
    public MVStore mvStore() {
        return mvStore;
    }

    @Bean("hg38variantAnnotator")
    @Override
    public VariantAnnotator variantAnnotator() {
        return super.buildVariantAnnotator();
    }

    @Bean("hg38variantFactory")
    @Override
    public VariantFactory variantFactory() {
        return super.buildVariantFactory();
    }

    @Bean("hg38variantDataService")
    @Override
    public VariantDataService variantDataService() {
        return super.buildVariantDataService();
    }

    @Bean("hg38genomeDataService")
    @Override
    public GenomeDataService genomeDataService() {
        return super.buildGenomeDataService();
    }

    //These require Spring to manage the caching and are called by buildVariantDataService
    @Bean("hg38genomeAnalysisService")
    @Override
    public GenomeAnalysisService genomeAnalysisService() {
        return buildGenomeAnalysisService();
    }

    @Bean("hg38allelePropertiesDao")
    @Override
    public AllelePropertiesDao allelePropertiesDao() {
        return new AllelePropertiesDaoMvStore(mvStore);
    }

    @Bean("hg38localFrequencyDao")
    @Override
    public FrequencyDao localFrequencyDao() {
        return new LocalFrequencyDao(localFrequencyTabixDataSource);
    }

    @Bean("hg38remmDao")
    @Override
    public RemmDao remmDao() {
        return new RemmDao(remmTabixDataSource);
    }

    @Bean("hg38caddDao")
    @Override
    public CaddDao caddDao() {
        return new CaddDao(caddIndelTabixDataSource, caddSnvTabixDataSource);
    }

    @Bean("hg38splicingDao")
    @Override
    public PathogenicityDao splicingDao() {
        ContigLengthDao contigLengthDao = new ContigLengthDao(splicingDataSource);
        SplicingPositionalWeightMatrixParser pwmParser = new DbSplicingPositionalWeightMatrixParser(splicingDataSource);
        SplicingTranscriptSource splicingTranscriptSource = new DbSplicingTranscriptSource(splicingDataSource, contigLengthDao.getContigLengths());

        GenomeCoordinatesFlipper flipper = new GenomeCoordinatesFlipper(contigLengthDao.getContigLengths());
        SplicingTranscriptLocator locator = new NaiveSplicingTranscriptLocator(pwmParser.getSplicingParameters(), flipper);
        SplicingInformationContentAnnotator annotator = new SplicingInformationContentAnnotator(pwmParser.getDonorMatrix(), pwmParser.getAcceptorMatrix(), pwmParser.getSplicingParameters());

        ScorerFactory scorerFactory = new RawScoringFactory(annotator);
        SplicingEvaluator evaluator = new SimpleSplicingEvaluator(locator, scorerFactory);

        return new SplicingDao(genomeSequenceAccessor, splicingTranscriptSource, evaluator);
    }

    @Bean("hg38testPathDao")
    @Override
    public PathogenicityDao testPathScoreDao() {
        return new TestPathogenicityScoreDao(testPathogenicitySource);
    }
}
