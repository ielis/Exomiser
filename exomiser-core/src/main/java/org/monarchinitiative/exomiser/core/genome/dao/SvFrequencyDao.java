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

package org.monarchinitiative.exomiser.core.genome.dao;

import org.monarchinitiative.exomiser.core.genome.ChromosomalRegionUtil;
import org.monarchinitiative.exomiser.core.model.ChromosomalRegion;
import org.monarchinitiative.exomiser.core.model.StructuralType;
import org.monarchinitiative.exomiser.core.model.Variant;
import org.monarchinitiative.exomiser.core.model.frequency.Frequency;
import org.monarchinitiative.exomiser.core.model.frequency.FrequencyData;
import org.monarchinitiative.exomiser.core.model.frequency.FrequencySource;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.sql.DataSource;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * @author Jules Jacobsen <j.jacobsen@qmul.ac.uk>
 */
public class SvFrequencyDao implements FrequencyDao {

    private final Logger logger = LoggerFactory.getLogger(SvFrequencyDao.class);

    private final DataSource svDataSource;

    public SvFrequencyDao(DataSource svDataSource) {
        this.svDataSource = svDataSource;
    }

    @Override
    public FrequencyData getFrequencyData(Variant variant) {
        if (variant.isStructuralVariant()) {
            int margin = ChromosomalRegionUtil.getBoundaryMargin(variant, 0.85);

            List<SvResult> results = runQuery(variant, margin);
//            logger.info("{}", variant);
//            logger.info("Searching for {}:{}-{}", variant.getChromosome(), variant.getStart() - margin, variant.getEnd() + margin);
//            results.forEach(svResult -> logger.info("{}", svResult));

            Map<Double, List<SvResult>> resultsByScore = results.stream()
                    .collect(Collectors.groupingBy(SvResult::getScore));
            List<SvResult> topMatches = resultsByScore.entrySet()
                    .stream()
                    .max(Map.Entry.comparingByKey())
                    .map(Map.Entry::getValue)
                    .orElse(List.of());

            if (topMatches.isEmpty()) {
                return FrequencyData.empty();
            }
            return mapToFrequencyData(topMatches);

        }
        return FrequencyData.empty();
    }

    private FrequencyData mapToFrequencyData(List<SvResult> topMatches) {
        SvResult first = topMatches.get(0);
        Frequency frequency = toFrequency(first);
        if (frequency.getFrequency() == 0) {
            // DGV has no frequency information, but does have an id
            return FrequencyData.of(first.id);
        }
        return FrequencyData.of(first.id, frequency);
    }

    private Frequency toFrequency(SvResult first) {
        FrequencySource frequencySource = frequencySource(first);
        return Frequency.of(frequencySource, first.af * 100);
    }

    private FrequencySource frequencySource(SvResult first) {
        switch (first.source) {
            case "GNOMAD_SV":
                return FrequencySource.GNOMAD_SV;
            case "DBVAR":
                return FrequencySource.DBVAR;
            case "DGV":
                return FrequencySource.DGV;
            case "GONL":
                return FrequencySource.GONL;
            case "DECIPHER":
                return FrequencySource.DECIPHER;
            default:
                return FrequencySource.UNKNOWN;
        }
    }

    private List<SvResult> runQuery(Variant variant, int margin) {
        String query = "SELECT * " +
                "FROM ( SELECT 'GNOMAD_SV' as SOURCE, " +
                "                CHR_ONE," +
                "                POS_ONE," +
                "                POS_TWO," +
                "                SV_LEN," +
                "                SV_TYPE," +
                "                ID," +
                "                AC," +
                "                AF" +
                "         FROM PBGA.GNOMAD_SV" +
                "         UNION ALL" +
                "         SELECT 'DBVAR' as SOURCE," +
                "                CHR_ONE," +
                "                POS_ONE," +
                "                POS_TWO," +
                "                SV_LEN," +
                "                SV_TYPE," +
                "                ID," +
                "                AC," +
                "                AF" +
                "         FROM PBGA.DBVAR_VARIANTS" +
                "         UNION ALL" +
                "         SELECT 'GONL' as SOURCE," +
                "                CHR_ONE," +
                "                POS_ONE," +
                "                POS_TWO," +
                "                SV_LEN," +
                "                SV_TYPE," +
                "                ID," +
                "                AN     as AC," +
                "                AF" +
                "         FROM PBGA.GONL" +
                "         UNION ALL" +
                "         SELECT 'DGV' as SOURCE," +
                "                CHR_ONE," +
                "                POS_ONE," +
                "                POS_TWO," +
                "                SV_LEN," +
                "                SV_TYPE," +
                "                ID," +
                "                AC," +
                "                0              as AF" +
                "         FROM PBGA.DGV_VARIANTS" +
                "         UNION ALL" +
                "         SELECT 'DECIPHER'    as SOURCE," +
                "                CHR_ONE," +
                "                POS_ONE," +
                "                POS_TWO," +
                "                SV_LEN," +
                "                SV_TYPE," +
                "                ID," +
                "                AC," +
                "                AF" +
                "         FROM PBGA.DECIPHER_CNV" +
                "     ) all_tables " +
                "WHERE CHR_ONE = ?" +
                "  and POS_ONE >= ? - ?" +
                "  and POS_ONE <= ? + ?" +
                "  and POS_TWO >= ? - ?" +
                "  and POS_TWO <= ? + ?" +
                "  and AC != -1;";

        try (Connection connection = svDataSource.getConnection();
             PreparedStatement ps = connection.prepareStatement(query)) {
            ps.setInt(1, variant.getChromosome());
            ps.setInt(2, variant.getStart());
            ps.setInt(3, margin);
            ps.setInt(4, variant.getStart());
            ps.setInt(5, margin);
            ps.setInt(6, variant.getEnd());
            ps.setInt(7, margin);
            ps.setInt(8, variant.getEnd());
            ps.setInt(9, margin);
            ResultSet rs = ps.executeQuery();

// consider also complex types where CHR_ONE != CHR_TWO - there are only about 600 in gnomad and gonl combined.
            return getSvResults(rs, variant);
        } catch (SQLException e) {
            logger.error("", e);
        }
        return List.of();
    }

    private List<SvResult> getSvResults(ResultSet rs, Variant variant) throws SQLException {

//            SOURCE	CHR_ONE	POS_ONE	POS_TWO	SV_LEN	SV_TYPE	ID	AC	AF
//            GNOMAD_SV	7	4972268	4973271	1003	DEL	gnomAD_v2_DEL_7_90956	94	0.004377
//            DBVAR	7	4972268	4973271	-1004	DEL	nssv15404135	94	0.004377
//            DGV_VARIANTS	7	4972233	4973326	1093	CNV	esv2659019	1151	0
//            DGV_VARIANTS	7	4972258	4973286	1028	CNV	esv3611958	2504	0
//            DECIPHER_CNV	7	4972251	4973293	1042	DEL	23156	94	0.049682875

        List<SvResult> results = new ArrayList<>();
        while (rs.next()) {
            String source = rs.getString("SOURCE");
            int chr = rs.getInt("CHR_ONE");
            int start = rs.getInt("POS_ONE");
            int end = rs.getInt("POS_TWO");
            int length = rs.getInt("SV_LEN");
            String svType = rs.getString("SV_TYPE");
            String id = rs.getString("ID");
            int ac = rs.getInt("AC");
            float af = rs.getFloat("AF");

            StructuralType structuralType = StructuralType.valueOf(svType);
            // there are cases such as INS_ME which won't match the database so we have to filter these here
            // consider also DEL/CNV_LOSS INS/CNV_GAIN/DUP/INS_ME and CNV
            if (structuralType.getBaseType() == variant.getStructuralType().getBaseType()) {
                SvResult svResult = new SvResult(chr, start, end, length, structuralType, source, id, ac, af);
                svResult.jaccard = ChromosomalRegionUtil.jaccard(variant, svResult);
//                System.out.println(svResult);
                results.add(svResult);
            }
        }
        return results;
    }

    private static class SvResult implements ChromosomalRegion {

        private final int chr;
        private final int start;
        private final int end;
        private final int length;
        private final StructuralType svType;
        private final String source;
        private final String id;
        private final int ac;
        private final float af;

        private double jaccard;

        private SvResult(int chr, int start, int end, int length, StructuralType svType, String source, String id, int ac, float af) {
            this.chr = chr;
            this.start = start;
            this.end = end;
            this.length = length;
            this.svType = svType;
            this.source = source;
            this.id = ".".equals(id) ? "" : id;
            this.ac = ac;
            this.af = af;
        }

        @Override
        public int getChromosome() {
            return chr;
        }

        @Override
        public int getStart() {
            return start;
        }

        @Override
        public int getEnd() {
            return end;
        }

        @Override
        public int getLength() {
            return length;
        }

        public double getScore() {
            // geometric mean of num alleles and similarity - try and get the best represented and most similar allele
            return Math.sqrt(ac * jaccard);
        }

        @Override
        public String toString() {
            return "SvResult{" +
                    "chr=" + chr +
                    ", start=" + start +
                    ", end=" + end +
                    ", length=" + length +
                    ", svType='" + svType + '\'' +
                    ", source='" + source + '\'' +
                    ", id='" + id + '\'' +
                    ", ac=" + ac +
                    ", af=" + af +
                    ", jaccard=" + jaccard +
                    ", score=" + getScore() +
                    '}';
        }
    }
}
