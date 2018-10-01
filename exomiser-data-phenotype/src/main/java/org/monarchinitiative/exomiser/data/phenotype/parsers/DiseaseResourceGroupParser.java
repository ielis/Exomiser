/*
 * The Exomiser - A tool to annotate and prioritize genomic variants
 *
 * Copyright (c) 2016-2017 Queen Mary University of London.
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

package org.monarchinitiative.exomiser.data.phenotype.parsers;

import org.monarchinitiative.exomiser.data.phenotype.resources.Resource;
import org.monarchinitiative.exomiser.data.phenotype.resources.ResourceGroup;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Component;

import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Jules Jacobsen <jules.jacobsen@sanger.ac.uk>
 */
@Component
public class DiseaseResourceGroupParser extends AbstractResourceGroupParser implements ResourceGroupParser {

    public static final String NAME = "DISEASETERMS";

    private static final Logger logger = LoggerFactory.getLogger(DiseaseResourceGroupParser.class);

    private Resource orphanet2GeneResource;
    private Resource disease2TermResource;

    @Override
    public void parseResources(ResourceGroup resourceGroup, Path inDir, Path outDir) {

        logger.info("Parsing {} resources. Writing out to: {}", resourceGroup.getName(), outDir);

        //Check everything is present before trying to parse them
        if (!requiredResourcesPresent(resourceGroup)) {
            logger.error("Not parsing {} ResourceGroup resources as not all required resources are present.", resourceGroup
                    .getName());
            return;
        }

        Map<String, String> disease2termMap = new HashMap<>();
        //first parseResource the mim2gene file
        Disease2TermParser disease2TermParser = new Disease2TermParser(disease2termMap);
        disease2TermParser.parseResource(disease2TermResource, inDir, outDir);

        Orphanet2GeneParser orphanet2GeneParser = new Orphanet2GeneParser(disease2termMap);
        orphanet2GeneParser.parseResource(orphanet2GeneResource, inDir, outDir);
    }

    @Override
    public boolean requiredResourcesPresent(ResourceGroup resourceGroup) {
        orphanet2GeneResource = resourceGroup.getResource(Orphanet2GeneParser.class);
        if (orphanet2GeneResource == null) {
            logResourceMissing(resourceGroup.getName(), Orphanet2GeneParser.class);
            return false;
        }

        disease2TermResource = resourceGroup.getResource(Disease2TermParser.class);
        if (disease2TermResource == null) {
            logResourceMissing(resourceGroup.getName(), Disease2TermParser.class);
            return false;
        }


        return true;
    }
}
