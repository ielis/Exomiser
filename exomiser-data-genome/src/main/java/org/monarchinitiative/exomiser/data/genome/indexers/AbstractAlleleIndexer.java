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

package org.monarchinitiative.exomiser.data.genome.indexers;

import org.monarchinitiative.exomiser.data.genome.model.Allele;
import org.monarchinitiative.exomiser.data.genome.model.AlleleResource;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.time.Duration;
import java.time.Instant;
import java.util.concurrent.atomic.AtomicLong;
import java.util.function.Consumer;

/**
 * @author Jules Jacobsen <j.jacobsen@qmul.ac.uk>
 */
public abstract class AbstractAlleleIndexer implements AlleleIndexer {

    private static final Logger logger = LoggerFactory.getLogger(AbstractAlleleIndexer.class);

    @Override
    public void index(AlleleResource alleleResource) {
        logger.info("Processing '{}' resource", alleleResource.getName());
        Instant startTime = Instant.now();
        AlleleLogger alleleLogger = new AlleleLogger(startTime);

        alleleResource.alleles()
                .peek(alleleLogger.logCount())
                .forEach(this::writeAllele);

        long seconds = Duration.between(startTime, Instant.now()).getSeconds();
        logger.info("Finished '{}' resource - processed {} alleles in {} sec. Total {} alleles written.",
                alleleResource.getName(),
                alleleLogger.count(),
                seconds,
                this.count());
    }

    protected abstract void writeAllele(Allele allele);

    public abstract long count();

    public abstract void close();

    private static class AlleleLogger {

        private final AtomicLong counter;
        private final Instant startTime;
        private Instant lastInstant;

        public AlleleLogger(Instant startTime) {
            this.counter = new AtomicLong();
            this.startTime = startTime;
            this.lastInstant = startTime;
        }

        public long count() {
            return counter.get();
        }

        public Consumer<Allele> logCount() {
            return allele -> {
                counter.incrementAndGet();
                int logInterval = 1000000;
                if (counter.get() % logInterval == 0) {
                    Instant now = Instant.now();
                    long totalSeconds = Duration.between(startTime, now).getSeconds();
                    long sinceLastCount = Duration.between(lastInstant, now).getSeconds();
                    long totalVar = counter.get();
                    logger.info("Processed {} variants total in {} sec - {} vps (last {} took {} sec - {} vps)", totalVar, totalSeconds, totalVar / totalSeconds, logInterval, sinceLastCount, logInterval / sinceLastCount);
                    logger.info("{}", allele);
                    lastInstant = now;
                }
            };
        }
    }

}
