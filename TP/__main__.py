from TP.loading import load_directory
from TP.kmers import stream_kmers, kmer2str

import os

def jaccard_from_sorted_lists(lstA, lstB):
    idxA = 0
    idxB = 0

    intersection = 0
    union = 0

    while idxA < len(lstA) and idxB < len(lstB):
        union += 1
        if lstA[idxA] == lstB[idxB]:
            intersection += 1
            idxA += 1
            idxB += 1
        elif lstA[idxA] < lstB[idxB]:
            idxA += 1
        else:
            idxB += 1

    union += len(lstA) - idxA
    union += len(lstB) - idxB

    return intersection / union



if __name__ == "__main__":
    print("Computation of Jaccard similarity between files")

    # Load all the files in a dictionary
    print("Loading files")
    files = load_directory("data")
    k = 21
    
    filenames = list(files.keys())
    
    # Create all the kmer lists (can be expensive in memory)
    print("Computing all kmer vectors")
    kmer_lists = {}
    for filename in filenames:
        kmer_lists[filename] = []
        # Enumerate all the sequences from a fasta
        for seq in files[filename]:
            kmer_lists[filename].extend(stream_kmers(seq, k))
        # Sort the kmer lists to speed up the comparison
        kmer_lists[filename].sort()

    print("Computing Jaccard similarity for all pairs of samples")
    for i in range(len(files)):
        for j in range(i+1, len(files)):
            jaccard = jaccard_from_sorted_lists(kmer_lists[filenames[i]], kmer_lists[filenames[j]])
            print(filenames[i], filenames[j], jaccard)

# def jaccard(fileA, fileB, k):
#     # j = 0
#     # --- To complete ---
#     index = stream_kmers(fileA, k)
#     intersect = 0
#     union = sum(index.values())
#     for kmer in enumerate_kmer(fileB, k):
#         if kmer in index and index[kmer] > 0:
#             intersect += 1
#             index[kmer] -= 1
#     else:
#         union += 1

#     return round(intersect/union, 3)
#     # return j


# if __name__ == "__main__":
#     print("Computation of Jaccard similarity between files")

#     # Load all the files in a dictionary
#     files = load_directory("data")
#     k = 21
            
#     print("Computing Jaccard similarity for all pairs of samples")
#     filenames = list(files.keys())
#     matrix = [list(range(0, len(filenames)+1))]
#     for i in range(1, len(filenames)+1):
#         matrix.append(list(range(i, len(filenames)+1+i)))
    
#     for i in range(len(files)):
#         matrix[0][i+1] = filenames[i]
#         matrix[i+1][0] = filenames[i]
#         for j in range(len(files)):
#             # --- Complete here ---
#             if i == j:
#                 matrix[i+1][j+1] = 0
#             else:
#                 J = jaccard(files[filenames[i]][0], files[filenames[j]][0], k)
#                 print(filenames[i], filenames[j], J)
#                 matrix[j+1][i+1] = J
    
#     print('Distance matrix: \n',matrix)

 
#     file_path = 'distance_matrix.txt'
#     if not os.path.exists(file_path):
#         with open(file_path, 'w') as file:
#             for row in matrix:
#                 file.write(str(row))
#                 file.write('\n')
#     else:
#         print(f"The file {file_path} already exists.")
