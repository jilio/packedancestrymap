package packedancestrymap

import (
	"bufio"
	"errors"
	"io"
	"math"
	"os"
	"regexp"
	"strconv"
	"strings"
	"sync"
)

// Snp from *.snp file
// 1st column is SNP name
// 2nd column is chromosome.  X chromosome is encoded as 23.
//     Also, Y is encoded as 24, mtDNA is encoded as 90, and XY is encoded as 91.
//     Note: SNPs with illegal chromosome values, such as 0, will be removed
// 3rd column is genetic position (in Morgans).  If unknown, ok to set to 0.0.
// 4th column is physical position (in bases)
// Optional 5th and 6th columns are reference and variant alleles.
//     For monomorphic SNPs, the variant allele can be encoded as X (unknown).
type Snp struct {
	Name        string
	Chromosome  uint8
	GeneticPos  float32
	PhysicalPos uint32
	Ref         uint8
	Alt         uint8
}

// Ind (individual) from *.ind file
// 1st column is sample ID.  Length is limited to 39 characters,
//     including the family name if that will be concatenated.
// 2nd column is gender (M or F).  If unknown, ok to set to U for Unknown.
// 3rd column is a label which might refer to Case or Control status,
//     or might be a population group label.
type Ind struct {
	SampleID string
	Gender   string
	Label    string
}

// ReadSnpFile - load *.snp file
func ReadSnpFile(path string) ([]Snp, error) {
	snps := []Snp{}
	file, err := os.Open(path)
	if err != nil {
		return snps, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		fields := strings.Fields(strings.Trim(scanner.Text(), " "))

		id := fields[0]
		chr, err := strconv.Atoi(fields[1])
		if err != nil {
			return snps, err
		}
		gpos, err := strconv.ParseFloat(fields[2], 32)
		if err != nil {
			return snps, err
		}
		ppos, err := strconv.ParseUint(fields[3], 10, 0)
		if err != nil {
			return snps, err
		}

		snps = append(snps, Snp{
			Name:        id,
			Chromosome:  uint8(chr),
			GeneticPos:  float32(gpos),
			PhysicalPos: uint32(ppos),
			Ref:         byte(fields[4][0]),
			Alt:         byte(fields[5][0]),
		})
	}
	if err := scanner.Err(); err != nil {
		return snps, err
	}

	return snps, nil
}

// ReadIndFile - load *.ind file
func ReadIndFile(path string) ([]Ind, error) {
	inds := []Ind{}

	file, err := os.Open(path)
	if err != nil {
		return inds, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		fields := strings.Fields(strings.Trim(scanner.Text(), " "))

		inds = append(inds, Ind{
			SampleID: fields[0],
			Gender:   fields[1],
			Label:    fields[2],
		})
	}
	if err := scanner.Err(); err != nil {
		return inds, err
	}

	return inds, nil
}

func ProcessGenoRows(genoPath, indPath, snpPath string, processFunc func(genoRow []byte, snp Snp, inds []Ind)) error {
	ok, err := Calcishash(genoPath, indPath, snpPath)
	if err != nil {
		return err
	}

	if !ok {
		return errors.New("Hash is not ok")
	}

	snps, err := ReadSnpFile(snpPath)
	if err != nil {
		return err
	}

	inds, err := ReadIndFile(indPath)
	if err != nil {
		return err
	}

	genoFile, err := os.Open(genoPath)
	if err != nil {
		return err
	}
	defer genoFile.Close()

	genoReader := bufio.NewReader(genoFile)
	chunkSize := int(math.Ceil(float64(len(inds)) / 4))

	// Why 48? See original EIG/src/mcio.c
	genoChunkSize := chunkSize
	if genoChunkSize < 48 {
		genoChunkSize = 48
	}

	// Skip header
	rchunk := make([]byte, genoChunkSize)
	genoReader.Read(rchunk)

	wg := &sync.WaitGroup{}
	for snpIndex := 0; snpIndex < len(snps); snpIndex++ {
		genoReader.Read(rchunk)

		genoRow := make([]byte, len(inds))
		for indIndex := 0; indIndex < len(inds); indIndex++ {
			byteOffset := int(math.Floor(float64(indIndex) / 4))
			bitOffset := uint8(indIndex%4) * 2
			b := rchunk[byteOffset]
			genotype := (b >> bitOffset) & 3
			genoRow[indIndex] = genotype
		}

		wg.Add(1)
		go func(genoRow []byte, snp Snp, inds []Ind, wg *sync.WaitGroup) {
			processFunc(genoRow, snp, inds)
			wg.Done()
		}(genoRow, snps[snpIndex], inds, wg)

	}
	wg.Wait()

	return nil
}

// HashIt takes a string as argument and calculates hash for it
func HashIt(str string) uint32 {
	var hash, length uint32 = 0, uint32(len(str))

	for j := uint32(0); j < length; j++ {
		hash *= 23
		hash += uint32(str[j])
	}

	return hash
}

// HashFileFirstColumn takes a file path string as argument
// and calculates hash for the file's first column (see *.snp or *.ind file)
// File must have last empty line
func HashFileFirstColumn(path string) (uint32, error) {
	file, err := os.Open(path)
	if err != nil {
		return 0, err
	}
	defer file.Close()

	var hash, thash uint32
	reader := bufio.NewReader(file)
	firstColumnRegexp, _ := regexp.Compile(`\S+`)
	for {
		str, err := reader.ReadString(10)
		if err != nil {
			if err == io.EOF {
				return hash, nil
			}
			return 0, err
		}

		thash = HashIt(firstColumnRegexp.FindString(str))
		hash *= 17
		hash ^= thash
	}
}

// Calcishash takes paths for *.geno, *.ind and *.snp files (eigenstrat combo)
// and calculate hashes on individuals and SNPs (to compare with file values).
func Calcishash(genoPath, indPath, snpPath string) (bool, error) {
	file, err := os.Open(genoPath)
	if err != nil {
		return false, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	ok := scanner.Scan()
	if !ok {
		return false, err
	}
	header := scanner.Text()

	hashInd, err := HashFileFirstColumn(indPath)
	if err != nil {
		return false, err
	}
	indHash := strconv.FormatInt(int64(hashInd), 16)

	hashSnp, err := HashFileFirstColumn(snpPath)
	if err != nil {
		return false, err
	}
	snpHash := strconv.FormatInt(int64(hashSnp), 16)

	indRegexp, err := regexp.Compile(indHash)
	if err != nil {
		return false, err
	}
	snpRegexp, err := regexp.Compile(snpHash)
	if err != nil {
		return false, err
	}
	hashOk := indRegexp.MatchString(header) && snpRegexp.MatchString(header) && (indPath != snpPath)

	return hashOk, nil
}
