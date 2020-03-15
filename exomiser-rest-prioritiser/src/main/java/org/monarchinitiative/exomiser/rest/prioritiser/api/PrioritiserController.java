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

package org.monarchinitiative.exomiser.rest.prioritiser.api;

import org.monarchinitiative.exomiser.core.genome.GenomeAnalysisService;
import org.monarchinitiative.exomiser.core.model.Gene;
import org.monarchinitiative.exomiser.core.model.GeneIdentifier;
import org.monarchinitiative.exomiser.core.prioritisers.HiPhiveOptions;
import org.monarchinitiative.exomiser.core.prioritisers.Prioritiser;
import org.monarchinitiative.exomiser.core.prioritisers.PriorityFactory;
import org.monarchinitiative.exomiser.core.prioritisers.PriorityResult;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.MediaType;
import org.springframework.web.bind.annotation.*;

import java.time.Duration;
import java.time.Instant;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static com.google.common.collect.ImmutableList.toImmutableList;

/**
 * @author Jules Jacobsen <jules.jacobsen@sanger.ac.uk>
 */
@RestController
public class PrioritiserController {

    private static final Logger logger = LoggerFactory.getLogger(PrioritiserController.class);

    private final PriorityFactory priorityFactory;
    private final Map<Integer, GeneIdentifier> geneIdentifiers;

    @Autowired
    public PrioritiserController(PriorityFactory priorityFactory, GenomeAnalysisService hg38GenomeAnalysisService) {
        this.priorityFactory = priorityFactory;
        Map<Integer, GeneIdentifier> map = new HashMap<>();
        for (GeneIdentifier geneIdentifier : hg38GenomeAnalysisService.getKnownGeneIdentifiers()) {
            // Don't add GeneIdentifiers without HGNC identifiers as these are superceeded by others with the same
            // entrez id which will create duplicate key errors and out of date gene symbols etc.
            if (geneIdentifier.hasEntrezId() && !geneIdentifier.getHgncId().isEmpty()) {
                GeneIdentifier previous = map.put(geneIdentifier.getEntrezIdAsInteger(), geneIdentifier);
                if (previous != null) {
                    logger.warn("Duplicate key added {} - was {}", geneIdentifier, previous);
                }
            }
        }
        this.geneIdentifiers = map;
        logger.info("Created GeneIdentifier cache with {} entries", geneIdentifiers.size());
    }

    @GetMapping(value = "about")
    public String about() {
        return "This service will return a collection of prioritiser results for any given set of:" +
                "\n\t - HPO identifiers e.g. HPO:00001" +
                "\n\t - Entrez gene identifiers e.g. 23364" +
                "\n\t - Specified prioritiser e.g. hiphive along with any prioritiser specific commands e.g. human,mouse,fish,ppi" +
                "\n\t - limit the number of genes returned e.g. 10";
    }

    @GetMapping(value = "", produces = MediaType.APPLICATION_JSON_UTF8_VALUE)
    public PrioritiserResultSet prioritise(@RequestParam(value = "phenotypes") Set<String> phenotypes,
                                           @RequestParam(value = "genes", required = false, defaultValue = "") Set<Integer> genesIds,
                                           @RequestParam(value = "prioritiser") String prioritiserName,
                                           @RequestParam(value = "prioritiser-params", required = false, defaultValue = "") String prioritiserParams,
                                           @RequestParam(value = "limit", required = false, defaultValue = "0") Integer limit
    ) {
        PrioritiserRequest prioritiserRequest = PrioritiserRequest.builder()
                .prioritiser(prioritiserName)
                .prioritiserParams(prioritiserParams)
                .genes(genesIds)
                .phenotypes(phenotypes)
                .limit(limit)
                .build();

        return prioritise(prioritiserRequest);
    }

    @PostMapping(value = "", consumes = MediaType.APPLICATION_JSON_UTF8_VALUE, produces = MediaType.APPLICATION_JSON_UTF8_VALUE)
    public PrioritiserResultSet prioritise(@RequestBody PrioritiserRequest prioritiserRequest) {
        logger.info("{}", prioritiserRequest);

        Instant start = Instant.now();

        Prioritiser<? extends PriorityResult> prioritiser = parsePrioritiser(prioritiserRequest.getPrioritiser(), prioritiserRequest
                .getPrioritiserParams());
        List<Gene> genes = makeGenesFromIdentifiers(prioritiserRequest.getGenes());

        List<PriorityResult> results = runLimitAndCollectResults(prioritiser, prioritiserRequest.getPhenotypes(), genes, prioritiserRequest
                .getLimit());

        Instant end = Instant.now();
        Duration duration = Duration.between(start, end);

        return new PrioritiserResultSet(prioritiserRequest, duration.toMillis(), results);
    }

    private Prioritiser<? extends PriorityResult> parsePrioritiser(String prioritiserName, String prioritiserParams) {
        switch (prioritiserName) {
            case "phenix":
                return priorityFactory.makePhenixPrioritiser();
            case "phive":
                return priorityFactory.makePhivePrioritiser();
            case "hiphive":
            default:
                HiPhiveOptions hiPhiveOptions = HiPhiveOptions.builder()
                        .runParams(prioritiserParams)
                        .build();
                return priorityFactory.makeHiPhivePrioritiser(hiPhiveOptions);
        }
    }

    private List<Gene> makeGenesFromIdentifiers(Collection<Integer> genesIds) {
        if (genesIds.isEmpty()) {
            logger.info("Gene identifiers not specified - will compare against all known genes.");
            //If not specified, we'll assume they want to use the whole genome. Should save people a lot of typing.
            //n.b. Gene is mutable so these can't be cached and returned.
            return allGenes();
        }
        // This is a hack - really the Prioritiser should only work on GeneIds, but currently this isn't possible as
        // OmimPrioritiser uses some properties of Gene
        return genesIds.stream()
                .map(id -> new Gene(geneIdentifiers.getOrDefault(id, unrecognisedGeneIdentifier(id))))
                .collect(toImmutableList());
    }

    private List<Gene> allGenes() {
        return geneIdentifiers.values().parallelStream()
                .map(Gene::new)
                .collect(toImmutableList());
    }

    private GeneIdentifier unrecognisedGeneIdentifier(Integer id) {
        return GeneIdentifier.builder().geneSymbol("GENE:" + id).build();
    }

    private <T extends PriorityResult> List<PriorityResult> runLimitAndCollectResults(Prioritiser<T> prioritiser, List<String> phenotypes, List<Gene> genes, int limit) {
        Set<Integer> wantedGeneIds = genes.stream().map(Gene::getEntrezGeneID).collect(Collectors.toSet());

        Stream<T> resultsStream = prioritiser.prioritise(phenotypes, genes)
                .filter(result -> wantedGeneIds.contains(result.getGeneId()))
                .sorted(Comparator.naturalOrder());

        logger.info("Finished {}", prioritiser.getPriorityType());
        if (limit == 0) {
            return resultsStream.collect(toImmutableList());
        }
        return resultsStream.limit(limit).collect(toImmutableList());
    }

}
