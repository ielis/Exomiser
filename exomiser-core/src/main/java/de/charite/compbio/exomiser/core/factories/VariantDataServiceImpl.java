/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package de.charite.compbio.exomiser.core.factories;

import de.charite.compbio.exomiser.core.dao.CaddDao;
import de.charite.compbio.exomiser.core.model.Variant;
import de.charite.compbio.exomiser.core.dao.FrequencyDao;
import de.charite.compbio.exomiser.core.dao.NcdsDao;
import de.charite.compbio.exomiser.core.dao.PathogenicityDao;
import de.charite.compbio.exomiser.core.dao.RegulatoryFeatureDao;
import de.charite.compbio.exomiser.core.model.frequency.FrequencyData;
import de.charite.compbio.exomiser.core.model.frequency.Frequency;
import de.charite.compbio.exomiser.core.model.frequency.FrequencySource;
import de.charite.compbio.exomiser.core.model.pathogenicity.PathogenicityData;
import de.charite.compbio.exomiser.core.model.pathogenicity.PathogenicityScore;
import de.charite.compbio.exomiser.core.model.pathogenicity.PathogenicitySource;
import de.charite.compbio.jannovar.annotation.VariantEffect;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import static java.util.stream.Collectors.toSet;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.context.annotation.Lazy;
import org.springframework.stereotype.Service;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 *
 * @author Jules Jacobsen <jules.jacobsen@sanger.ac.uk>
 */
@Service
public class VariantDataServiceImpl implements VariantDataService {

    private static final Logger logger = LoggerFactory.getLogger(VariantDataServiceImpl.class);

    @Autowired
    private FrequencyDao frequencyDao;
    @Autowired
    private PathogenicityDao pathogenicityDao;
    @Lazy
    @Autowired
    private CaddDao caddDao;
    @Lazy
    @Autowired
    private NcdsDao ncdsDao;
    @Autowired
    private RegulatoryFeatureDao regulatoryFeatureDao;
    
    private static final PathogenicityData NO_PATH_DATA = new PathogenicityData();

    @Override
    public FrequencyData getVariantFrequencyData(Variant variant, Set<FrequencySource> frequencySources) {
        FrequencyData allFrequencyData = frequencyDao.getFrequencyData(variant);
        return frequencyDataFromSpecifiedSources(allFrequencyData, frequencySources);
    }

    protected FrequencyData frequencyDataFromSpecifiedSources(FrequencyData allFrequencyData, Set<FrequencySource> frequencySources) {
        Set<Frequency> wanted = allFrequencyData.getKnownFrequencies().stream()
                .filter(frequency -> frequencySources.contains(frequency.getSource()))
                .collect(toSet());
        return new FrequencyData(allFrequencyData.getRsId(), wanted);
    }

    @Override
    public PathogenicityData getVariantPathogenicityData(Variant variant, Set<PathogenicitySource> pathogenicitySources) {
        //OK, this is a bit stupid, but if no sources are defined we're not going to bother checking for data
        if (pathogenicitySources.isEmpty()) {
            return NO_PATH_DATA;
        }
        //TODO: ideally we'd have some sort of compact, high-performance document store for this sort of data rather than several different datasources query and ship.
        List<PathogenicityScore> allPathScores = new ArrayList<>();
        //Polyphen, Mutataion Taster and SIFT are all trained on missense variants - this is what is contained in the original variant table, but we shouldn't know that.
        if (variant.getVariantEffect() == VariantEffect.MISSENSE_VARIANT) {
            PathogenicityData missenseScores = pathogenicityDao.getPathogenicityData(variant);
            allPathScores.addAll(missenseScores.getPredictedPathogenicityScores());
        } else if (pathogenicitySources.contains(PathogenicitySource.NCDS)) {
        //NCDS is trained on all the non-coding bits of the genome, this outperforms CADD for non-coding variants
            PathogenicityData nonCodingScore = ncdsDao.getPathogenicityData(variant);
            allPathScores.addAll(nonCodingScore.getPredictedPathogenicityScores());
        }
        
        //CADD does all of it although is not as good as NCDS for the non-coding regions.
        if (pathogenicitySources.contains(PathogenicitySource.CADD)) {
            PathogenicityData caddScore = caddDao.getPathogenicityData(variant);
            allPathScores.addAll(caddScore.getPredictedPathogenicityScores());
        }

        return pathDataFromSpecifiedDataSources(allPathScores, pathogenicitySources);
    }

    protected PathogenicityData pathDataFromSpecifiedDataSources(List<PathogenicityScore> allPathScores, Set<PathogenicitySource> pathogenicitySources) {
        Set<PathogenicityScore> wanted = allPathScores.stream()
                .filter(pathogenicity -> pathogenicitySources.contains(pathogenicity.getSource()))
                .collect(toSet());
        return new PathogenicityData(wanted);
    }
    
    @Override
    public VariantEffect getVariantRegulatoryFeatureData(Variant variant) {
        if (variant.getVariantEffect() == VariantEffect.INTERGENIC_VARIANT || variant.getVariantEffect() == VariantEffect.UPSTREAM_GENE_VARIANT) {
            return regulatoryFeatureDao.getRegulatoryFeatureData(variant);
        }
        return variant.getVariantEffect();
    }

}
