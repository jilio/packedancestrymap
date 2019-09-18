package packedancestrymap

import "testing"

func TestPackage(t *testing.T) {
	datasets := []string{
		"test_data/Anc_EA2",
		"test_data/Olalde_et_al_genotypes",
		"test_data/Roopkund_1240k",
	}

	t.Run("Hash check", func(t *testing.T) {
		for _, dataset := range datasets {
			geno := dataset + ".geno"
			ind := dataset + ".ind"
			snp := dataset + ".snp"

			got, err := Calcishash(geno, ind, snp)
			want := true
			if err != nil {
				t.Error(err)
			}

			if got != want {
				t.Errorf("\n%s:\n  got %t\n  want %t", dataset, got, want)
			}
		}
	})
}
