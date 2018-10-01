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

package org.monarchinitiative.exomiser.core.genome;

import java.util.Objects;

/**
 * genome reference assembly version - hg19/hg38.
 */
public enum GenomeAssembly {

    HG19("hg19"), HG38("hg38");

    private final String value;

    GenomeAssembly(String value) {
        this.value = value;
    }

    public static GenomeAssembly defaultBuild() {
        return GenomeAssembly.HG19;
    }

    public static GenomeAssembly fromValue(String value) {
        Objects.requireNonNull(value, "Genome build cannot be null");
        switch (value.toLowerCase()) {
            case "hg19":
            case "hg37":
            case "grch37":
                return HG19;
            case "hg38":
            case "grch38":
                return HG38;
            default:
                throw new InvalidGenomeAssemblyException(String.format("'%s' is not a valid/supported genome assembly.", value));
        }
    }

    @Override
    public String toString() {
        return value;
    }

    public static class InvalidGenomeAssemblyException extends RuntimeException {

        public InvalidGenomeAssemblyException() {
        }

        public InvalidGenomeAssemblyException(String message) {
            super(message);
        }

        public InvalidGenomeAssemblyException(String message, Throwable cause) {
            super(message, cause);
        }

        public InvalidGenomeAssemblyException(Throwable cause) {
            super(cause);
        }

        public InvalidGenomeAssemblyException(String message, Throwable cause, boolean enableSuppression, boolean writableStackTrace) {
            super(message, cause, enableSuppression, writableStackTrace);
        }
    }
}