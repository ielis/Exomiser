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

package org.monarchinitiative.exomiser.core.model.frequency;

import com.fasterxml.jackson.annotation.JsonIgnore;
import com.google.common.collect.Maps;

import java.util.*;

/**
 * Frequency data for the variant from the Thousand Genomes, the Exome Server
 * Project and Broad ExAC datasets.
 *
 * Note that the frequency data are expressed as percentages.
 *
 * @author Jules Jacobsen <jules.jacobsen@sanger.ac.uk>
 */
public class FrequencyData {

    private static final FrequencyData EMPTY_DATA = new FrequencyData(RsId.empty(), Collections.emptyMap());

    private static final float VERY_RARE_SCORE = 1f;
    private static final float NOT_RARE_SCORE = 0f;

    private final RsId rsId;
    private final Map<FrequencySource, Frequency> knownFrequencies;

    public static FrequencyData of(RsId rsId, Collection<Frequency> frequencies) {
        return validate(rsId, frequencies);
    }

    public static FrequencyData of(RsId rsId, Frequency frequency) {
        return validate(rsId, Collections.singletonList(frequency));
    }

    public static FrequencyData of(RsId rsId, Frequency... frequency) {
        return validate(rsId, Arrays.asList(frequency));
    }

    public static FrequencyData of(Frequency frequency) {
        return validate(RsId.empty(), Collections.singletonList(frequency));
    }

    public static FrequencyData of(Frequency... frequency) {
        return validate(RsId.empty(), Arrays.asList(frequency));
    }

    public static FrequencyData empty() {
        return EMPTY_DATA;
    }

    private static FrequencyData validate(RsId rsId, Collection<Frequency> frequencies) {
        Objects.requireNonNull(rsId, "RsId cannot be null");
        Objects.requireNonNull(frequencies, "frequencies cannot be null");

        if (rsId.isEmpty() && frequencies.isEmpty()) {
            return FrequencyData.empty();
        }
        Map<FrequencySource, Frequency> frequencySourceMap = new EnumMap<>(FrequencySource.class);
        for (Frequency frequency : frequencies) {
            frequencySourceMap.put(frequency.getSource(), frequency);
        }
        return new FrequencyData(rsId, frequencySourceMap);
    }

    private FrequencyData(RsId rsId, Map<FrequencySource, Frequency> knownFrequencies) {
        this.rsId = rsId;
        this.knownFrequencies = Maps.immutableEnumMap(knownFrequencies);
    }

    //RSID ought to belong to the Variant, not the frequencyData, but its here for convenience
    public RsId getRsId() {
        return rsId;
    }

    public Frequency getFrequencyForSource(FrequencySource source) {
        return knownFrequencies.get(source);
    }

    /**
     * @return true if this variant is at all represented in dbSNP or ESP data,
     * regardless of frequency. That is, if the variant has an RS id in dbSNP or
     * any frequency data at all, return true, otherwise false.
     */
    @JsonIgnore
    public boolean isRepresentedInDatabase() {
        return hasDbSnpRsID() || hasKnownFrequency();
    }

    public boolean hasDbSnpData() {
        return knownFrequencies.containsKey(FrequencySource.THOUSAND_GENOMES);
    }

    public boolean hasDbSnpRsID() {
        return !rsId.isEmpty();
    }

    public boolean hasEspData() {
        for (FrequencySource dataSource : knownFrequencies.keySet()) {
            switch (dataSource) {
                case ESP_AFRICAN_AMERICAN:
                case ESP_EUROPEAN_AMERICAN:
                case ESP_ALL:
                    return true;
                default:
            }
        }
        return false;
    }
    
    public boolean hasExacData() {
        for (FrequencySource dataSource : knownFrequencies.keySet()) {
            switch (dataSource) {
                case EXAC_AFRICAN_INC_AFRICAN_AMERICAN:
                case EXAC_AMERICAN:
                case EXAC_EAST_ASIAN:
                case EXAC_FINNISH:
                case EXAC_NON_FINNISH_EUROPEAN:
                case EXAC_OTHER:
                case EXAC_SOUTH_ASIAN:
                    return true;
                default:
            }
        }
        return false;
    }

    public boolean hasKnownFrequency() {
        return !knownFrequencies.isEmpty();
    }

    /**
     * This function tests whether or not this {@code FrequencyData} object contains a {@code Frequency} object which has
     * a frequency greater than the maximum frequency provided. This method does not check any ranges so it is advised
     * that the user checks the frequency type in advance of calling this method. By default exomiser expresses the
     * frequencies as a <b>percentage</b> value.
     *
     * @param maxFreq the maximum frequency threshold against which the {@code Frequency} objects are tested
     * @return true if the object contains a {@code Frequency} over the provided percentage value, otherwise returns false.
     * @since 10.1.0
     */
    public boolean hasFrequencyOverPercentageValue(float maxFreq) {
        for (Frequency frequency : knownFrequencies.values()) {
            if (frequency.isOverThreshold(maxFreq)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Returns a list of {@code Frequency} objects. If there is no known frequency data then an empty list will be returned.
     * This method will return a mutable copy of the underlying data.
     *
     * @return a mutable copy of the {@code Frequency} data
     */
    public List<Frequency> getKnownFrequencies() {
        return new ArrayList<>(knownFrequencies.values());
    }

    /**
     * Returns a the maximum frequency - if there are no known frequencies/ no
     * frequency data it will return 0.
     *
     * @return
     */
    @JsonIgnore
    public float getMaxFreq() {
        return (float) knownFrequencies.values().stream().mapToDouble(Frequency::getFrequency).max().orElse(0);
    }

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 29 * hash + Objects.hashCode(this.rsId);
        hash = 29 * hash + Objects.hashCode(this.knownFrequencies);
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final FrequencyData other = (FrequencyData) obj;
        if (!Objects.equals(this.rsId, other.rsId)) {
            return false;
        }
        return Objects.equals(this.knownFrequencies, other.knownFrequencies);
    }

    @Override
    public String toString() {
        return "FrequencyData{" + "rsId=" + rsId + ", knownFrequencies=" + knownFrequencies.values() + '}';
    }

    /**
     * @return returns a numerical value that is closer to one, the rarer
     * the variant is. If a variant is not entered in any of the data
     * sources, it returns one (highest score). Otherwise, it identifies the
     * maximum MAF in any of the databases, and returns a score that depends on
     * the MAF. Note that the frequency is expressed as a percentage.
     */
    public float getScore() {

        float max = getMaxFreq();

        if (max <= 0) {
            return VERY_RARE_SCORE;
        } else if (max > 2) {
            return NOT_RARE_SCORE;
        } else {
            return 1.13533f - (0.13533f * (float) Math.exp(max));
        }
    }

}
