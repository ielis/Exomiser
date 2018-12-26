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

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.monarchinitiative.exomiser.core.writers;

import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.SerializationFeature;
import de.charite.compbio.jannovar.mendel.ModeOfInheritance;
import org.monarchinitiative.exomiser.core.analysis.Analysis;
import org.monarchinitiative.exomiser.core.analysis.AnalysisResults;
import org.monarchinitiative.exomiser.core.model.Gene;
import org.monarchinitiative.exomiser.core.prioritisers.HiPhivePriorityResult;
import org.monarchinitiative.exomiser.core.prioritisers.PriorityType;
import org.monarchinitiative.exomiser.core.writers.phenogrid.PhenoGrid;
import org.monarchinitiative.exomiser.core.writers.phenogrid.PhenoGridAdaptor;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Jules Jacobsen <jules.jacobsen@sanger.ac.uk>
 */
public class PhenogridWriter implements ResultsWriter {

    private static final Logger logger = LoggerFactory.getLogger(PhenogridWriter.class);    
    private static final OutputFormat OUTPUT_FORMAT = OutputFormat.PHENOGRID;
    
    @Override
    public void writeFile(ModeOfInheritance modeOfInheritance, Analysis analysis, AnalysisResults analysisResults, OutputSettings settings) {
        String outFileName = ResultsWriterUtils.makeOutputFilename(analysis.getVcfPath(), settings.getOutputPrefix(), OUTPUT_FORMAT, ModeOfInheritance.ANY);
        Path outFile = Paths.get(outFileName);
        try (BufferedWriter writer = Files.newBufferedWriter(outFile, Charset.defaultCharset())) {
            //HiPhive phenotype scores are not dependant on the mode of inheritance.
            writer.write(writeString(ModeOfInheritance.ANY, analysis, analysisResults, settings));

        } catch (IOException ex) {
            logger.error("Unable to write results to file {}", outFileName, ex);
        }
        logger.debug("{} results written to file {}", OUTPUT_FORMAT, outFileName);

    }

    @Override
    public String writeString(ModeOfInheritance modeOfInheritance, Analysis analysis, AnalysisResults analysisResults, OutputSettings settings) {
        List<Gene> passedGenes = ResultsWriterUtils.getMaxPassedGenes(analysisResults.getGenes(), settings.getNumberOfGenesToShow());
        List<HiPhivePriorityResult> hiPhiveResults = new ArrayList<>();
        for (Gene gene : passedGenes) {
            if (gene.getPriorityResults().containsKey(PriorityType.HIPHIVE_PRIORITY)) {
               HiPhivePriorityResult hiPhiveResult = (HiPhivePriorityResult) gene.getPriorityResult(PriorityType.HIPHIVE_PRIORITY);
               hiPhiveResults.add(hiPhiveResult);
            }
        }
        PhenoGridAdaptor phenoGridAdaptor = new PhenoGridAdaptor();
        PhenoGrid phenogrid = phenoGridAdaptor.makePhenoGridFromHiPhiveResults("hiPhive specified phenotypes", hiPhiveResults) ;
        return writePhenoGridAsJson(phenogrid);
    }

    private String writePhenoGridAsJson(PhenoGrid phenogrid) {
        ObjectMapper mapper = new ObjectMapper();
        mapper.configure(SerializationFeature.INDENT_OUTPUT, false);
        mapper.configure(SerializationFeature.WRITE_ENUMS_USING_TO_STRING, true);
        try {
            return mapper.writeValueAsString(phenogrid);
        } catch (JsonProcessingException ex) {
            logger.error("Error mapping PhenoGrid to json", ex);
        }
        return "Error mapping PhenoGrid to json";
    }

}
