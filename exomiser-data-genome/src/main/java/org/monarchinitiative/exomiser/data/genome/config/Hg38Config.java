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

package org.monarchinitiative.exomiser.data.genome.config;

import com.google.common.collect.ImmutableMap;
import org.monarchinitiative.exomiser.data.genome.model.AlleleResource;
import org.monarchinitiative.exomiser.data.genome.model.parsers.DbNsfpColumnIndex;
import org.monarchinitiative.exomiser.data.genome.model.resource.*;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.core.env.Environment;

import java.net.URL;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Map;

/**
 * @author Jules Jacobsen <j.jacobsen@qmul.ac.uk>
 */
@Configuration
public class Hg38Config {

    @Autowired
    public Environment environment;

    @Bean
    public Path hg38GenomePath() {
        return Paths.get(environment.getProperty("hg38.genome-dir"));
    }

    @Bean
    public Map<String, AlleleResource> hg38AlleleResources() {
        ImmutableMap.Builder<String, AlleleResource> alleleResources = new ImmutableMap.Builder<>();

        alleleResources.put("gnomad-genome", gnomadGenomeAlleleResource());
        alleleResources.put("gnomad-exome", gnomadExomeAlleleResource());
        // TOPMed removed as this is now part of dbSNP
        alleleResources.put("dbsnp", dbSnpAlleleResource());
        alleleResources.put("uk10k", uk10kAlleleResource());
        alleleResources.put("exac", exacAlleleResource());
        alleleResources.put("esp", espAlleleResource());
        alleleResources.put("dbnsfp", dbnsfpAlleleResource());
        alleleResources.put("clinvar", clinVarAlleleResource());

        return alleleResources.build();
    }

    public AlleleResource dbSnpAlleleResource() {
        String namespacePrefix = "hg38.dbsnp";
        AlleleResourceProperties resourceProperties = getAlleleResourceProperties(namespacePrefix);
        Path resourcePath = resourceProperties.getAlleleResourcePath();
        URL resourceUrl = resourceProperties.getAlleleResourceUrl();
        return new DbSnpAlleleResource(namespacePrefix, resourceUrl, resourcePath);
    }

    public AlleleResource clinVarAlleleResource() {
        String namespacePrefix = "hg38.clinvar";
        AlleleResourceProperties resourceProperties = getAlleleResourceProperties(namespacePrefix);
        Path resourcePath = resourceProperties.getAlleleResourcePath();
        URL resourceUrl = resourceProperties.getAlleleResourceUrl();
        return new ClinVarAlleleResource(namespacePrefix, resourceUrl, resourcePath);
    }

    public AlleleResource espAlleleResource() {
        String namespacePrefix = "hg38.esp";
        AlleleResourceProperties resourceProperties = getAlleleResourceProperties(namespacePrefix);
        Path resourcePath = resourceProperties.getAlleleResourcePath();
        URL resourceUrl = resourceProperties.getAlleleResourceUrl();
        return new EspHg38AlleleResource(namespacePrefix, resourceUrl, resourcePath);
    }

    public AlleleResource exacAlleleResource() {
        String namespacePrefix = "hg38.exac";
        AlleleResourceProperties resourceProperties = getAlleleResourceProperties(namespacePrefix);
        Path resourcePath = resourceProperties.getAlleleResourcePath();
        URL resourceUrl = resourceProperties.getAlleleResourceUrl();
        return new ExacExomeAlleleResource(namespacePrefix, resourceUrl, resourcePath);
    }

    public AlleleResource dbnsfpAlleleResource() {
        String namespacePrefix = "hg38.dbnsfp";
        AlleleResourceProperties resourceProperties = getAlleleResourceProperties(namespacePrefix);
        Path resourcePath = resourceProperties.getAlleleResourcePath();
        URL resourceUrl = resourceProperties.getAlleleResourceUrl();
        if (resourcePath.toString().contains("dbNSFP4.")) {
            return new DbNsfp4AlleleResource(namespacePrefix, resourceUrl, resourcePath, DbNsfpColumnIndex.HG38);
        }
        return new DbNsfp3AlleleResource(namespacePrefix, resourceUrl, resourcePath, DbNsfpColumnIndex.HG38);
    }

    public AlleleResource topmedAlleleResource() {
        String namespacePrefix = "hg38.topmed";
        AlleleResourceProperties resourceProperties = getAlleleResourceProperties(namespacePrefix);
        Path resourcePath = resourceProperties.getAlleleResourcePath();
        URL resourceUrl = resourceProperties.getAlleleResourceUrl();
        return new TopMedAlleleResource(namespacePrefix, resourceUrl, resourcePath);
    }

    public AlleleResource uk10kAlleleResource() {
        String namespacePrefix = "hg38.uk10k";
        AlleleResourceProperties resourceProperties = getAlleleResourceProperties(namespacePrefix);
        Path resourcePath = resourceProperties.getAlleleResourcePath();
        URL resourceUrl = resourceProperties.getAlleleResourceUrl();
        return new Uk10kAlleleResource(namespacePrefix, resourceUrl, resourcePath);
    }

    public AlleleResource gnomadGenomeAlleleResource() {
        String namespacePrefix = "hg38.gnomad-genome";
        AlleleResourceProperties resourceProperties = getAlleleResourceProperties(namespacePrefix);
        Path resourcePath = resourceProperties.getAlleleResourcePath();
        URL resourceUrl = resourceProperties.getAlleleResourceUrl();
        return new GnomadGenomeAlleleResource(namespacePrefix, resourceUrl, resourcePath);
    }

    public AlleleResource gnomadExomeAlleleResource() {
        String namespacePrefix = "hg38.gnomad-exome";
        AlleleResourceProperties resourceProperties = getAlleleResourceProperties(namespacePrefix);
        Path resourcePath = resourceProperties.getAlleleResourcePath();
        URL resourceUrl = resourceProperties.getAlleleResourceUrl();
        return new GnomadExomeAlleleResource(namespacePrefix, resourceUrl, resourcePath);
    }

    private AlleleResourceProperties getAlleleResourceProperties(String namespacePrefix) {
        String fileName = environment.getProperty(namespacePrefix + ".file-name");
        Path fileDir = Paths.get(environment.getProperty(namespacePrefix + ".file-dir"));
        String fileUrl = environment.getProperty(namespacePrefix + ".file-url");
        return new AlleleResourceProperties(fileName, fileDir, fileUrl);
    }
}
